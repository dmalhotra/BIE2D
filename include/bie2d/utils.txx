#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

namespace sctl {

  void commandline_option_start(int argc, char** argv, const char* help_text, const Comm& comm) {
    if (comm.Rank()) return;
    char help[] = "--help";
    for (int i = 0; i < argc; i++) {
      if (!strcmp(argv[i], help)) {
        if (help_text != NULL) std::cout << help_text << std::endl;
        std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
      }
    }
  }

  const char* commandline_option(int argc, char** argv, const char* opt, const char* default_val, bool required, bool implicit, const char* desc, const Comm& comm){
    char help[] = "--help";
    for (int i = 0; i < argc; i++) {
      if (!strcmp(argv[i], help)) {
        if (!comm.Rank()) {
          std::cout << "  " << std::left << std::setw(14) << opt;
          std::cout << "  " << (desc ? desc : "");
          std::cout << " (" << default_val << ")";
          std::cout << std::endl;
        }
        return default_val;
      }
    }

    for (int i = 0; i < argc; i++) {
      if (!strcmp(argv[i], opt)) {
        if (implicit) {
          return "true";
        } else {
          return argv[(i+1)%argc];
        }
      }
    }
    if (required) {
      if (!comm.Rank()) std::cout << "Missing required option\n" << "    " << opt << "  " << (desc?desc:"") << "\n\n";
      if (!comm.Rank()) std::cout << "To see usage options\n" << "    "<<argv[0]<<" --help\n\n";
      exit(0);
    }
    return default_val;
  }

  void commandline_option_end(int argc, char** argv) {
    char help[] = "--help";
    for (int i = 0; i < argc; i++) {
      if (!strcmp(argv[i], help)) {
        exit(0);
      }
    }
  }

  bool strtob(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
  }

}
