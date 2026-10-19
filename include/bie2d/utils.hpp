// Copyright (C) 2023-2026 The Simons Foundation, Inc. and the BIE2D authors (see NOTICE)
// Authors: Dhairya Malhotra, Dan Fortunato, Mariana Martinez Aguilar
// SPDX-License-Identifier: Apache-2.0

#ifndef BIE2D_UTILS
#define BIE2D_UTILS

#include <sctl.hpp>

namespace sctl {

    void commandline_option_start(int argc, char** argv, const char* help_text = nullptr, const Comm& comm = Comm::Self());

    const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, bool implicit, const char* err_msg, const Comm& comm = Comm::Self());

    void commandline_option_end(int argc, char** argv);

    bool strtob(std::string str);

}

#include <bie2d/utils.txx>

#endif

