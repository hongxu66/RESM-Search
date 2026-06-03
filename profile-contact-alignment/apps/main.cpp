/*************************************************************************
    > File Name: main.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun May 15 16:20:23 2022
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

//    CmapSearch(OPT);
    CmapAli(OPT);
    //if(OPT.msa){map2seqswithmsa(OPT);}
    //else if(OPT.gen) mulmap2genome(OPT);
    //else mulmap2mulseqs(OPT);

    return 0;
}
