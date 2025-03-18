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

