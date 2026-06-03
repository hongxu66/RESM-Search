/*************************************************************************
    > File Name: buildmap.h
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun May 15 16:20:00 2022
 ************************************************************************/

#ifndef __BUILDMAP_H
#define __BUILDMAP_H

#include "global_variable.h"
#include "load_data.h"
#include "msa.h"
#include "buildmap.h"
//#include "msa.h"


class Map{
    public:
        mtx_double mtx;
        vec_int vec;
        vec_int vec_div;
        mtx_int vec_i; 
        mtx_double prf; 
        vec_char aa; 
        vec_char ss;
        vec_int M2N;
        string id;
        vec_string os_map;
        vec_string os_prf;
//        string sequence;
        Map(){};
        Map(string ID):id(ID){};
        //vec_int load_data_file(string file,string &seq, int sep_cutoff, vec_bool &range,float prob, string type);
//        vec_int load_data(string file, int sep_cutoff);
        vec_int load_data(string file,string &seq, int sep_cutoff);
        vec_int load_data_from_file(string file,string &seq, int sep_cutoff, float prob, string type);
        //vec_int seq2data(string seq,int sep_cutoff, vec_bool &range, float prob,string type);
        vec_int seq2data(string seq,int sep_cutoff, float prob,string type);
        //vec_int load_data(string file, int sep_cutoff, vec_bool &range);

        vec_int buildmapsMSA(string &seq,string &seqName, vec_string msanames,vec_string msaalns, int sep_cutoff);
        vec_int buildmapsMSA(string &seq,string &seqName,string &seqss,vec_string msanames,vec_string msalans, int sep_cutoff);
        vec_int buildmaps(string &seq,string &ss, int sep_cutoff,float prob);
        vec_int buildmapswithMSA(string &seq,string &seqName,string &seqss,vec_string msanames,vec_string msalans, int sep_cutoff);
        int printout(string outfile);
        ~Map(){}
};

vec_string seq2prf(string seq, float prob=0.85);
vec_string seq2prf(string seq, float prob, string ss);
vec_string seq2map(string seq);
vec_string ss2map(string ss);

void alignseqs(string seq1,string seq2, vector<vector<int>> maps);


Matrix4d PAM(int num);

//void Options::get_opt(vec_string opt);


void mod_gap(vec_double &gap, vec_char &ss, double &gap_ss_w);

sequence * readsinglefasta(const std::string filename);

#endif
