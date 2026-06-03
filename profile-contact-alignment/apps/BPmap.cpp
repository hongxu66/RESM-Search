/*************************************************************************
    BPmap: align two precomputed contact map files.
 ************************************************************************/

#include "align.h"

using namespace std;

namespace {
string basename_of(const string &path) {
    size_t pos = path.find_last_of("/\\");
    return pos == string::npos ? path : path.substr(pos + 1);
}

void ensure_sequence(string &seq, const Map &map) {
    if(seq.empty()) {
        seq.assign(map.mtx.size(), 'N');
    }
}

int align_contact_maps(Options &OPT) {
    string seqa;
    string seqb;

    Map A(basename_of(OPT.filea));
    vec_int m2n_a = A.load_data(OPT.filea, seqa, OPT.sep_cutoff);
    ensure_sequence(seqa, A);

    int size_a = A.mtx.size();
    if(size_a == 0) {
        cout << "ERROR: empty contact map: " << OPT.filea << endl;
        return 1;
    }

    vec_double gap_a(size_a, OPT.gap_open);
    if(OPT.use_gap_ss == true) {
        mod_gap(gap_a, A.ss, OPT.gap_ss_w);
    }

    Map B(basename_of(OPT.fileb));
    vec_int m2n_b = B.load_data(OPT.fileb, seqb, OPT.sep_cutoff);
    ensure_sequence(seqb, B);

    int size_b = B.mtx.size();
    if(size_b == 0) {
        cout << "ERROR: empty contact map: " << OPT.fileb << endl;
        return 1;
    }

    vec_double gap_b(size_b, OPT.gap_open);
    if(OPT.use_gap_ss == true) {
        mod_gap(gap_b, B.ss, OPT.gap_ss_w);
    }

    mtx_double P_SCO;
    if(OPT.use_prf == true) {
        ini_prf_SCO(P_SCO, OPT.prf_w, B, A);
    }

    Alignment Ali(B, A, OPT, seqb, seqa);
    vec_int a2b_max = Ali.mapalignment(gap_b, gap_a, P_SCO);

    return a2b_max.empty() ? 1 : 0;
}
}

int main(int argc, const char *argv[]) {
    clock_t start_time = clock();
    set_resm_run_start(start_time);

    Options OPT;
    vec_string opt;
    for(int a = 1; a < argc; a++) {
        string arg = argv[a];
        opt.push_back(arg.c_str());
    }
    OPT.get_opt(opt);

    int status = align_contact_maps(OPT);

    return status;
}
