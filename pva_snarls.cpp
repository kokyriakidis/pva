// pva_snarls.cpp — `pva snarls` subcommand
// Wraps: vg snarls <gfa> > <out.snarls>

#include "pva.h"
#include "pva_utils.h"

#include <iostream>
#include <string>

namespace pva {

int snarls(const std::string& gfa,
           const std::string& out_snarls,
           const std::string& vg_bin)
{
    std::string bin = pva::utils::find_bin(vg_bin, "vg", "vg");
    int ret = pva::utils::run_to_file({bin, "snarls", gfa}, out_snarls);
    if (ret != 0)
        std::cerr << "pva snarls: vg snarls failed (exit " << ret << ")\n";
    return ret;
}

static void snarls_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " snarls [options]\n"
        "\n"
        "Options:\n"
        "  --gfa <path>      Input GFA file\n"
        "  --out <path>      Output snarls file\n"
        "  --vg-bin <path>   Path to vg binary (auto-detected if omitted)\n"
        "\n"
        "Example:\n"
        "  pva snarls --gfa region.gfa --out region.snarls\n";
}

int snarls_main(int argc, char* argv[]) {
    std::string gfa, out_snarls, vg_bin;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva snarls: " << a << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if      (a == "--gfa")               gfa        = next();
        else if (a == "--out")               out_snarls = next();
        else if (a == "--vg-bin")            vg_bin     = next();
        else if (a == "--help" || a == "-h") { snarls_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva snarls: unknown option " << a << "\n";
            snarls_usage(argv[0]); return 1;
        }
    }

    if (gfa.empty() || out_snarls.empty()) { snarls_usage(argv[0]); return 1; }
    return pva::snarls(gfa, out_snarls, vg_bin);
}

} // namespace pva
