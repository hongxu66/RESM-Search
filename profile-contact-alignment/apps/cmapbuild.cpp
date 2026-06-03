/*************************************************************************
    > File Name: cmapbuild.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun Oct 23 11:27:16 2022
 ************************************************************************/

//#include "global_variable.h"
//#include "buildmap.h"
#include "align.h"
//#include "msa.h"
//#include "testmap.h"

using namespace std;

void print_help(){
    cout << "${cmap_bin}/cmap_build\n";
    cout << "                      \n";
    cout << "Options: -dbn dbnfile \n";
    cout << "         -sto stockholm\n";
    cout << "         -a2m a2mfile \n";
    cout << "         -o output.CMP\n";
    cout << "example: ${cmap_bin}/cmap_build -dbn dbnfile -o output.CMP\n";
}



int main(int argc, char * argv[])
{
    clock_t start_time = clock();
    set_resm_run_start(start_time);

    Options OPT;
    vec_string opt;
    for(int a=1;a<argc;a++){string arg=argv[a];opt.push_back(arg.c_str());}

    string outputfile{""};
    string dbnfile{""};
    string a2mfile{""};
    string stofile{""};

    if(count(opt.begin(),opt.end(),"-dbn") && count(opt.begin(),opt.end(),"-sto") && count(opt.begin(),opt.end(),"-a2m")) {print_help();exit(1);}

    int c;
    static struct option long_options[] = {
        {"dbn",required_argument,NULL,'d'},
        {"a2m",no_argument,      NULL,'a'},
        {"sto",no_argument,      NULL,'s'},
        {"out",no_argument,      NULL,'o'}
    };

    while(1){
        int opt_index=0;
        c = getopt_long(argc,argv,"d:a:s:o:",long_options,&opt_index);
        if(-1 == c){
            break;
        }
        switch(c){
            case 'd':
                dbnfile = optarg;
//                cout << optarg << endl;
//                cout << "dbnfile:" << dbnfile << endl;
                break;
            case 'a':
                a2mfile = optarg;
//                cout << "a2mfile:" << a2mfile << endl;
                break;
            case 's':
                stofile = optarg;
//                cout << "stofile:" << stofile << endl;
                break;
            case 'o':
                outputfile = optarg;
//                cout << "outputfile:" << outputfile << endl;
                break;
            default:
                print_help();
                exit(1);
//                break;
        }
    }





//    string outname;
//
////    cout << "argc:" << argc;
//    if(argc <= 3) {print_help();exit(1);}
//    for(int a=1;a<argc;a++){string arg=argv[a];opt.push_back(arg.c_str());}
//
//
//    if(count(arg.begin(),arg.end(),"-out")){
//        outputfile = arg[find(arg.begin(),arg.end(),"-out") - arg.begin()];
//    }
//    else{
//        outputfile = "out.CMP";
//    }
//
//    if(count(arg.begin(),arg.end(),"-dbn")){
//        dbnfile = arg[find(arg.begin(),arg.end(),"-dbn") - arg.begin()];
//    }
//    else if(count(arg.begin(),arg.end(),"-a2m")){
//        a2mfile = arg[find(arg.begin(),arg.end(),"-a2m") - arg.begin()];
//    }
//    else if(count(arg.begin(),arg.end(),"-sto")){
//        stofile = arg[find(arg.begin(),arg.end(),"-sto") - arg.begin()];
//    }
//    else {print_help();exit(1);}
//
//    int flag = 0;
//    for(int i=0;i<opt.size();i++){
//        string arg=opt[i];
////        cout << arg << endl;
//        if(arg == "-dbnfile") {i++;flag++;}
//        else if(arg == "-msafile") {i++;flag++;}
//        else if(arg == "-stofile") {i++;flag=3;}
//        else if(arg=="-o") {outname = opt[i+1];i++;}
//        else continue;
//    }
//    OPT.get_opt(opt);
//
    Map A;
    seqdbn *queryseqss;
    if(dbnfile.length()) {
        cout << "dbnfile:" << dbnfile;
        queryseqss = readdbnfile(dbnfile);
    }
    else {print_help();exit(1);}

    if(outputfile.length()==0){outputfile="TMP.cmp";}

//    vec_string aNames, aSeqs, aSS;
    if(a2mfile.length()){
        vec_string msaNames, msaAlns;
        readfastafile(msaNames,msaAlns,a2mfile);
        if(msaAlns.size()==0){cout<<"MSAfile is empty!!!"<<endl;exit(1);}
        vec_int m2n_a = A.buildmapswithMSA(queryseqss->seq,queryseqss->id,queryseqss->ss,msaNames,msaAlns,OPT.sep_cutoff);
    }
    else if(stofile.length()){
        vec_string msaNames, msaAlns;
        string alnRNAtype,alnRNAss;
        int alnNum = readMSAsto(stofile,msaNames,msaAlns,alnRNAtype,alnRNAss);
        cout << "ALNNUM:" << alnNum << endl;
        if(alnNum==0) {cout << "MSAfile is empty!!!"<<endl;exit(1);}
        vec_int m2n_a = A.buildmapswithMSA(queryseqss->seq,queryseqss->id,queryseqss->ss,msaNames,msaAlns,OPT.sep_cutoff);
    }
    else{
        A.buildmaps(queryseqss->seq,queryseqss->ss,OPT.sep_cutoff,OPT.prob);
    }

////    os_map = ss2map(queryseqss->ss);
////    os_prf = seq2prf(queryseqss->seq,OPT.prob,queryseqss->ss);
//    if(flag==1) {}
//    else if(flag==2) {;}
//    else if(flag==3) {;}
    A.printout(outputfile);

//    cout << "Query:" << queryseqss->id << endl;
//    seq = queryseqss->seq;


    return 0;
}
