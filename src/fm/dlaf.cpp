//
// Distributed Linear Algebra with Future (DLAF)
//
// Copyright (c) 2018-2023, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#include <dlaf/init.h>
#include <pika/execution.hpp>
#include <pika/init.hpp>
#include <pika/program_options.hpp>
#include <pika/runtime.hpp>

static bool dlaf_init_ = false;

extern "C" void dlaf_init() {
  if (!dlaf_init_) {
    int argc = 1;
    const char *const argv[] = {"cp2k", nullptr};

    pika::program_options::options_description desc("cp2k");
    desc.add(dlaf::getOptionsDescription());

    /* pika initialization */
    pika::init_params p;
    p.rp_callback = dlaf::initResourcePartitionerHandler;
    p.desc_cmdline = desc;
    pika::start(nullptr, argc, argv, p);

    /* DLA-Future initialization */
    dlaf::initialize(argc, argv);
    dlaf_init_ = true;
    pika::suspend();
  }
}
