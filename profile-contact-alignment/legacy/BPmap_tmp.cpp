/*************************************************************************
    > File Name: RESM_Search_align (BPmap entry: BPmap_tmp.cpp)
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun Oct 23 11:27:26 2022
 ************************************************************************/

#include <time.h>
//#include "global_variable.h"
//#include "buildmap.h"
#include "BPmap.h"
//#include "msa.h"
//#include "testmap.h"

using namespace std;

int main(int argc, const char * argv[])
{
    clock_t start_time,end_time;
    start_time = clock();


    Options OPT;

    vec_string opt;for (int a = 1; a < argc; a++){string arg = argv[a];opt.push_back(arg.c_str());}
    OPT.get_opt(opt);

//    CmapAli(OPT);
    CmapAlign(OPT);
    //if(OPT.msa){map2seqswithmsa(OPT);}
    //else if(OPT.gen) mulmap2genome(OPT);
    //else mulmap2mulseqs(OPT);

    end_time = clock();
    double run_time = double (end_time - start_time) / CLOCKS_PER_SEC;
    cout << "RUN TIME:" << run_time << endl;

    return 0;
}
