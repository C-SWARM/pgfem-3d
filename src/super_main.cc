/*** This is the PGFem3D main function. It has a branch between the
     macro-micro or single scale execution models. */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <map>

extern int single_scale_main(int, char *[]);
extern int multi_scale_main(int, char *[]);

static constexpr const char USAGE[] =
 R"(PGFem3D.

     Usage:
       PGFem3D (-SS | -MS) {options}

     Options:
       -SS, --SS   Single scale simulation
       -MS, --SS   Multi-scale simulation
       -h, --help  Show this screen
)";

/// Map from valid options to handlers (help just dumps usage).
static std::map<std::string, decltype(&single_scale_main)> frameworks = {
  { "-SS",    single_scale_main },
  { "--SS",   single_scale_main },
  { "-MS",    multi_scale_main },
  { "--MS",   multi_scale_main },
  { "-h",     [](int, char*[]) { std::cout << USAGE; return 0; } },
  { "--help", [](int, char*[]) { std::cout << USAGE; return 0; } }
};

/// Dispatch to the correct framework.
int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "ERROR: must provide at least one option!\n\n" << USAGE;
  }
  else if (auto main = frameworks[argv[1]]) {
    return main(argc, argv);
  }
  else {
    std::cerr << "ERROR: unrecognized framework " << argv[1] << "\n\n" << USAGE;
  }
  return 1;
}

