/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_gpu_kernel.h"
#include "dbm_multiply_opencl.cl.h"

typedef struct {
  int max_m, max_n, avg_m, avg_n, avg_k, changes;
} dbm_multiply_gpu_launch_info_t;

static void dbm_multiply_gpu_launch_info(dbm_multiply_gpu_launch_info_t *info,
                                         const dbm_task_t *tasks, int ntasks) {
  int i = 1;
  info->max_m = tasks[0].m;
  info->avg_m = tasks[0].m;
  info->max_n = tasks[0].n;
  info->avg_n = tasks[0].n;
  info->avg_k = tasks[0].k;
  for (info->changes = 0; i < ntasks; ++i) {
    const int m = tasks[i].m, n = tasks[i].n, k = tasks[i].k;
    info->max_m = imax(info->max_m, m);
    info->max_n = imax(info->max_n, n);
    if (info->avg_m != m || info->avg_n != n || info->avg_k != k) {
      info->avg_m = (info->avg_m + m) / 2;
      info->avg_n = (info->avg_n + n) / 2;
      info->avg_k = (info->avg_k + k) / 2;
      ++info->changes;
    }
  }
}

void dbm_multiply_gpu_launch_kernel(const offloadStream_t stream, double alpha,
                                    int ntasks, const dbm_task_t *tasks_host,
                                    const dbm_task_t *tasks,
                                    const double *pack_a_data,
                                    const double *pack_b_data,
                                    double *shard_c_data) {
  /* creating/calling kernel must be consistent across threads */
  static cl_kernel kernel_global = NULL;
  static LIBXSMM_TLS cl_kernel kernel = NULL;
  static int ndims = 1, clinear = 0;
  static size_t wgsize[] = {0, 0, 0};
  const libxsmm_timer_tickint start = libxsmm_timer_tick();
  const c_dbcsr_acc_opencl_config_t *const config = &c_dbcsr_acc_opencl_config;
  const int verbosity = config->verbosity;
  int result = EXIT_SUCCESS;
  cl_event event, *const perf_event =
                      ((0 <= verbosity && 2 >= verbosity) ? NULL : &event);
  const c_dbcsr_acc_opencl_stream_t *const str = ACC_OPENCL_STREAM(stream);
  size_t work_size[] = {1, 1, 1}, ibatch = 0;
  size_t iadata = 0, ibdata = 0, icdata = 0;
  const size_t work_tasks = ntasks;
  dbm_multiply_gpu_launch_info_t info = {0};
  c_dbcsr_acc_opencl_info_memptr_t adata, bdata, cdata, batch;
  assert(NULL != pack_a_data && NULL != pack_b_data && NULL != shard_c_data);
  assert(NULL != str && NULL != str->queue);
  assert(0 < ntasks && NULL != tasks);
#if defined(OPENCL_DBM_SOURCE_MULTIPLY)
  if (NULL == kernel_global) { /* initial check if kernel is present */
    ACC_OPENCL_ACQUIRE(config->lock_main);
    if (NULL == kernel_global) {
      char params[ACC_OPENCL_BUFFERSIZE] =
          "-cl-fast-relaxed-math -cl-denorms-are-zero";
      const char *const gen_env = getenv("DBM_MULTIPLY_GEN");
      const char *const lin_env = getenv("DBM_MULTIPLY_LIN");
      const char *const bn_env = getenv("DBM_MULTIPLY_BN");
      const char *const sm_env = getenv("DBM_MULTIPLY_SM");
      const char *const wg_env = getenv("DBM_MULTIPLY_WG");
      const char *const lu_env = getenv("DBM_MULTIPLY_LU");
      const char *const xf_env = getenv("DBM_MULTIPLY_XF");
      const c_dbcsr_acc_opencl_device_t *const devinfo = &config->device;
      int sm = (NULL == sm_env ? 0 /*default*/ : atoi(sm_env));
      const int bn0 = (0 == devinfo->nv ? (0 == devinfo->amd ? 4 : 8) : 2);
      const int bn1 = ((0 == sm && 0 == clinear) ? bn0 : (bn0 * 2));
      int bn = LIBXSMM_CLMP(NULL == bn_env ? bn1 : atoi(bn_env), 1, 32);
      int lu = LIBXSMM_CLMP(NULL == lu_env ? 0 : atoi(lu_env), -2, 1);
      int gen = ((NULL == bn_env && NULL == sm_env && NULL == wg_env &&
                  NULL == lu_env && NULL == lin_env)
                     ? (NULL == gen_env ? 1 /*default*/ : atoi(gen_env))
                     : 0);
      const int gpu = (CL_DEVICE_TYPE_GPU == devinfo->type);
      const int xf = (NULL == xf_env ? -1 /*default*/ : atoi(xf_env));
      const char *extensions[] = {NULL, NULL}, *flags = NULL;
      size_t nextensions = sizeof(extensions) / sizeof(*extensions);
      const size_t wgsize0 = devinfo->wgsize[0], wgsize1 = devinfo->wgsize[1];
      size_t wgsize2 = devinfo->wgsize[2];
      size_t offset =
          ((0 == config->debug && 0 == config->dump) ? strlen(params) : 0);
      offset += (size_t)c_dbcsr_acc_opencl_flags_atomics(
          devinfo, c_dbcsr_acc_opencl_atomic_fp_64, extensions, &nextensions,
          params + offset, sizeof(params) - offset);
      if (2 <= gen || (0 != gen && 0 != wgsize2 /*subgroups*/ &&
                       2 <= *devinfo->std_level && NULL != extensions[1] &&
                       NULL != strstr(extensions[1], "cl_ext_float_atomics"))) {
        offset +=
            (size_t)LIBXSMM_SNPRINTF(params + offset, sizeof(params) - offset,
                                     " -DDBM_MULTIPLY_OPENCL_GEN");
        wgsize[1] = wgsize[2] = 1;
        wgsize[0] = 16;
        lu = bn = 0;
        ndims = 3;
      } else {
        wgsize[0] = (NULL == wg_env ? (unsigned long int)LIBXSMM_ABS(sm)
                                    : strtoul(wg_env, NULL, 10));
        if (0 != wgsize2 && 0 < wgsize[0]) { /* subgroups */
          if (LIBXSMM_DELTA(wgsize[0], wgsize1) <=
              LIBXSMM_DELTA(wgsize[0], wgsize2)) { /* select SG-size */
            wgsize2 = wgsize1;
          }
          wgsize[0] = LIBXSMM_UP(wgsize[0], wgsize2);
        } else {
          wgsize[0] = LIBXSMM_UP(wgsize[0], wgsize1);
          wgsize2 = 0;
        }
        wgsize[0] = LIBXSMM_CLMP(wgsize[0], 0, wgsize0);
        sm = ((0 != sm && 0 != wgsize[0])
                  ? (LIBXSMM_ISPOT(bn * sizeof(double)) + 1)
                  : 0);
        clinear = (NULL == lin_env ? 0 /*default*/ : atoi(lin_env));
        offset += (size_t)LIBXSMM_SNPRINTF(
            params + offset, sizeof(params) - offset,
            " %s %s -DBN=%i -DSM=%i -DLU=%i -DWG=%i -DSG=%i",
            0 != gpu ? "-DGPU" : "", 0 == clinear ? "" : "-DCLINEAR", bn, sm,
            lu, (int)wgsize[0], (int)wgsize2);
        gen = 0;
      }
      if (0 != devinfo->intel && 0 < xf) {
        flags = "-cl-intel-256-GRF-per-thread";
      }
      result |= (sizeof(params) > offset ? EXIT_SUCCESS : EXIT_FAILURE);
      result |= c_dbcsr_acc_opencl_kernel(
          0 /*source_is_file*/, OPENCL_DBM_SOURCE_MULTIPLY, "dbm_multiply",
          params, flags, NULL /*try*/, NULL /*try_ok*/, extensions, nextensions,
          &kernel_global);
      if (2 <= verbosity || 0 > verbosity) {
        if (EXIT_SUCCESS == result) {
          const double ds = libxsmm_timer_duration(start, libxsmm_timer_tick());
          fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel gpu=%i", gpu);
          if (0 != gen) { /* generated kernel */
            fprintf(stderr, " gen=%i", gen);
          }
          if (0 != clinear) {
            fprintf(stderr, " lin=%i", clinear);
          }
          if (0 != bn) {
            fprintf(stderr, " bn=%i", bn);
          }
          if (0 != sm) {
            fprintf(stderr, " sm=%i", sm);
          }
          if (0 != wgsize[0]) {
            fprintf(stderr, " wg=%i", (int)wgsize[0]);
          }
          if (0 != wgsize2) {
            fprintf(stderr, " sg=%i", (int)wgsize2);
          }
          if (0 != lu) {
            fprintf(stderr, " lu=%i", lu);
          }
          fprintf(stderr, " ms=%.1f\n", 1E3 * ds);
        } else {
          fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel failed to generate\n");
        }
      }
    }
    kernel = clCloneKernel(kernel_global, &result); /* always clone */
    ACC_OPENCL_RELEASE(config->lock_main);
  } else if (NULL == kernel) {
    kernel = clCloneKernel(kernel_global, &result);
  }
#else
#error "OpenCL kernel code not found!"
#endif
  result |= c_dbcsr_acc_opencl_info_devptr_lock(&adata, NULL /*lock*/,
                                                pack_a_data, 1 /*esize*/,
                                                NULL /*amount*/, &iadata);
  result |= c_dbcsr_acc_opencl_info_devptr_lock(&bdata, NULL /*lock*/,
                                                pack_b_data, 1 /*esize*/,
                                                NULL /*amount*/, &ibdata);
  result |= c_dbcsr_acc_opencl_info_devptr_lock(&cdata, NULL /*lock*/,
                                                shard_c_data, 1 /*esize*/,
                                                NULL /*amount*/, &icdata);
  result |= c_dbcsr_acc_opencl_info_devptr_lock(
      &batch, NULL /*lock*/, tasks /*batch*/, sizeof(dbm_task_t), &work_tasks,
      &ibatch);
  assert(0 == iadata && 0 == ibdata && 0 == icdata);
  result |= clSetKernelArg(kernel, 0, sizeof(cl_double), &alpha);
  result |= clSetKernelArg(kernel, 1, sizeof(cl_int), &ibatch);
  if (1 < ndims) { /* DBM_MULTIPLY_GEN */
    const cl_uint zero = 0;
    assert(0 != wgsize[1] && 0 != wgsize[1] && 0 != wgsize[2]);
    work_size[0] = 16;
    assert(1 == work_size[1]);
    work_size[2] = work_tasks;
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 2, batch.memory);
    result |= clSetKernelArg(kernel, 3, sizeof(cl_uint), &zero /*shape*/);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 4, adata.memory);
    result |= clSetKernelArg(kernel, 5, sizeof(cl_uint), &zero /*A_shape0*/);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, bdata.memory);
    result |= clSetKernelArg(kernel, 7, sizeof(cl_uint), &zero /*B_shape0*/);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 8, cdata.memory);
    result |= clSetKernelArg(kernel, 9, sizeof(cl_uint), &zero /*C_shape0*/);
  } else {
    size_t size = work_tasks;
    dbm_multiply_gpu_launch_info(&info, tasks_host, ntasks);
    size *= (0 == clinear ? info.max_m : info.max_n);
    /* fixup to be a multiple of the WG-size */
    work_size[0] = (0 < wgsize[0] ? LIBXSMM_UP(size, wgsize[0]) : size);
    result |= clSetKernelArg(kernel, 2, sizeof(cl_int), &ntasks);
    result |= clSetKernelArg(kernel, 3, sizeof(cl_int), &size);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 4, batch.memory);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 5, adata.memory);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, bdata.memory);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 7, cdata.memory);
  }
  result |= clEnqueueNDRangeKernel(
      str->queue, kernel, ndims, NULL, work_size, 0 < wgsize[0] ? wgsize : NULL,
      0 /*num_wait*/, NULL /*wait_list*/, perf_event);
  if (NULL != perf_event && EXIT_SUCCESS == result &&
      EXIT_SUCCESS == clWaitForEvents(1, perf_event)) {
    const double dhost = libxsmm_timer_duration(start, libxsmm_timer_tick());
    cl_ulong begin = 0, end = 0;
    if (EXIT_SUCCESS ==
            clGetEventProfilingInfo(*perf_event, CL_PROFILING_COMMAND_START,
                                    sizeof(cl_ulong), &begin, NULL) &&
        EXIT_SUCCESS == clGetEventProfilingInfo(*perf_event,
                                                CL_PROFILING_COMMAND_END,
                                                sizeof(cl_ulong), &end, NULL)) {
      const double dkrnl = 1E-9 * LIBXSMM_DELTA(begin, end);
      const double dtotl = 1E+3 * LIBXSMM_MAX(dkrnl, dhost);
      int pure;
      if (1 < ndims) { /* DBM_MULTIPLY_GEN */
        dbm_multiply_gpu_launch_info(&info, tasks_host, ntasks);
      }
      pure = (100 * (ntasks - info.changes) + ntasks - 1) / ntasks;
      fprintf(stderr,
              "INFO ACC/LIBDBM: DBM-kernel mnk=%ix%ix%i pure=%i%% "
              "ntasks=%i kernel_ms=%.2g total_ms=%.2g gflops=%.1f\n",
              info.avg_m, info.avg_n, info.avg_k, pure, ntasks, dkrnl, dtotl,
              2E-6 * info.avg_m * info.avg_n * info.avg_k * ntasks / dtotl);
    }
  }
  OFFLOAD_CHECK(result);
}

#endif // defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

// EOF
