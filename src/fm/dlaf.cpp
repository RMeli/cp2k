//
// Distributed Linear Algebra with Future (DLAF)
//
// Copyright (c) 2018-2023, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#include <dlaf_c/init.h>

extern "C" void dlaf_init() {
  const char* pika_argv[] = {"cp2k", "--pika:print-bind"};
  const char* dlaf_argv[] = {"cp2k"};
  dlaf_initialize(2, pika_argv, 1, dlaf_argv);
}
