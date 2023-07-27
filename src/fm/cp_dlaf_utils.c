#if defined(__DLAF)

#include <dlaf_c/init.h>

void dlaf_init() {
  const char *pika_argv[] = {"cp2k", "--pika:print-bind"};
  const char *dlaf_argv[] = {"cp2k"};
  dlaf_initialize(2, pika_argv, 1, dlaf_argv);
}

#endif // __DLAF
