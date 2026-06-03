/*************************************************************************
    > File Name: cmapsearch.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun Oct 23 11:27:43 2022
 ************************************************************************/

//#include "global_variable.h"
//#include "buildmap.h"
#include "align.h"
//#include "msa.h"
//#include "testmap.h"

using namespace std;

void print_help(){
    cout << "${cmap_bin}/cmap_search\n";
    cout << "                      \n";
    cout << "Options: -dbn dbnfile \n";
    cout << "         -sto stockholm\n";
    cout << "         -a2m a2mfile \n";
    cout << "         -cmap cmapfile \n";
    cout << "         -genome fastafile \n";
    cout << "         -database fastafile \n";
    cout << "example: ${cmap_bin}/cmap_search {-dbn dbnfile} {-dbn dbnfile -sto stofile} {-dbn dbnfile -a2m a2mfile} {-camp cmapfile} {-genome fastafile} {-database fastafile}\n";
}


int main(int argc, char * argv[])
{
    clock_t start_time = clock();
    set_resm_run_start(start_time);


    Options OPT;

    vec_string opt;for (int a = 1; a < argc; a++){string arg = argv[a];opt.push_back(arg.c_str());}
//    OPT.get_opt(opt);

    string outputfile{""};
    string dbnfile{""};
    string a2mfile{""};
    string stofile{""};
    string cmapfile("");
    string genomefile{""};
    string databasefile{""};

    if(count(opt.begin(),opt.end(),"-dbn") && count(opt.begin(),opt.end(),"-sto") && count(opt.begin(),opt.end(),"-a2m") && count(opt.begin(),opt.end(),"-cmap")) {print_help();exit(1);}
    if(count(opt.begin(),opt.end(),"-database") && count(opt.begin(),opt.end(),"-genome") ) {print_help();exit(1);}


    int c;
    static struct option long_options[] = {
        {"dbn",required_argument,NULL,'d'},
        {"a2m",no_argument,      NULL,'a'},
        {"sto",no_argument,      NULL,'s'},
        {"cmap",no_argument,      NULL,'c'},
        {"genome",no_argument,    NULL,'g'},
        {"database",no_argument,  NULL,'b'}
//        {"out",no_argument,      NULL,'o'}
    };

    while(1){
        int opt_index=0;
        c = getopt_long(argc,argv,"d:a:s:c:g:b:",long_options,&opt_index);
        if(-1 == c){
            break;
        }
        switch(c){
            case 'd':
                dbnfile = optarg;
//                cout << optarg << endl;
                cout << "dbnfile:" << dbnfile << endl;
                break;
            case 'a':
                a2mfile = optarg;
//                cout << "a2mfile:" << a2mfile << endl;
                break;
            case 's':
                stofile = optarg;
//                cout << "stofile:" << stofile << endl;
                break;
//            case 'o':
//                outputfile = optarg;
//                cout << "outputfile:" << outputfile << endl;
                break;
            case 'c':
                cmapfile = optarg;
                break;
            case 'b':
                databasefile = optarg;
                break;
            case 'g':
                genomefile = optarg;
                break;
            default:
                print_help();
                exit(1);
//                break;
        }
    }


    if(cmapfile.length()){
        cout << "Underdeveloped !!!";
        exit(1);
    }
    else if(dbnfile.length()){
        OPT.filea = dbnfile;
        if(stofile.length()) {OPT.msa=true; OPT.msafile=stofile;}
        else if(a2mfile.length()) {OPT.msa=true; OPT.msafile=a2mfile;/*cout<<"Underdeveloped"<<endl;exit(1);*/}
        if(genomefile.length()) {OPT.fileb = genomefile; CmapSearchGenome(OPT);}
        else if(databasefile.length()) {OPT.fileb = databasefile;CmapSearchDatabase(OPT);}
        else {print_help();exit(1);}
    }
    else {print_help();exit(1);}


//    CmapSearch(OPT);
    //if(OPT.msa){map2seqswithmsa(OPT);}
    //else if(OPT.gen) mulmap2genome(OPT);
    //else mulmap2mulseqs(OPT);

    return 0;
}
