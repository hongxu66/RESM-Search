/*************************************************************************
    > File Name: RESM_Search_alignNCBI.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun Oct 23 11:27:26 2022
 ************************************************************************/

//#include "global_variable.h"
//#include "buildmap.h"
#include "align.h"
//#include "msa.h"
//#include "testmap.h"

using namespace std;

int main(int argc, const char * argv[])
{
    clock_t start_time = clock();
    set_resm_run_start(start_time);


    Options OPT;

    vec_string opt;for (int a = 1; a < argc; a++){string arg = argv[a];opt.push_back(arg.c_str());}
    OPT.get_opt(opt);

//    CmapAli(OPT);
//    CmapAlignNCBI(OPT);
//    CmapAlignNCBI1201(OPT);
    CmapAlignNCBI1202(OPT);
    //if(OPT.msa){map2seqswithmsa(OPT);}
    //else if(OPT.gen) mulmap2genome(OPT);
    //else mulmap2mulseqs(OPT);

    return 0;
}
