/*************************************************************************
    > File Name: load_data.h
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Thu Oct 20 11:06:36 2022
 ************************************************************************/
#ifndef __LOAD_DATA_H
#define __LOAD_DATA_H

#include "global_variable.h"

struct seqdbn{
    string id;
    string seq;
    string ss;
};
struct sequence{
    string id;
    string seq;
};


class Options{
    public:
        // for cmapbuild
        string dbnfile;
        string stofile;
        // for cmapsearch
        string filea,fileb;
//        string seqa,seqb;
//        vec_bool range_a,range_b;
        bool use_gap_ss = true;
        double gap_ss_w = 1;
//        bool use_prf = true;
        bool use_prf = false;
        double prf_w = 0;
        float prob = 0.6;
        double gap_open = -1;
//        double gap_open = -1;
//        double gap_ext = -0.00;
        double gap_ext = -0.01;
        int sep_cutoff = 1;
        int iter = 20;
        bool silent = false;
        bool msa = false;
        bool gen = false;
//        int lengthC = 10000;
        int lengthC = 1000;
        int lengthE = 1000;
        int overlap = 200;
        int numThreads = 8;
        string msafile;
        string msafileb;

        void get_opt(vec_string opt);
};


seqdbn * readdbnfile(const std::string filename);

string readctfile(const std::string filename);

sequence * readsinglefasta(const std::string filename);

void readdbnfile(vec_string &seqnames, vec_string &seqs, vec_string &seqss,const std::string filename);
void readfastafile(vec_string & refnames, vec_string & refseqs, const std::string filename);
void readfastafile(Options OPT, vec_string & refnames, vec_string & refseqs, const std::string filename);



#endif
