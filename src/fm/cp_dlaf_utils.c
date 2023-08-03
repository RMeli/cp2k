/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2023 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__DLAF)

#include <dlaf_c/init.h>

void dlaf_init() {
  const char *pika_argv[] = {"cp2k", "--pika:bind=none", "--pika:threads=12"};
  const char *dlaf_argv[] = {"cp2k"};
  dlaf_initialize(3, pika_argv, 1, dlaf_argv);
}

#endif // __DLAF
