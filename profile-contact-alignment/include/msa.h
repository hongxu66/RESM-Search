/*************************************************************************
    > File Name: msa.h
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun May 15 16:19:22 2022
 ************************************************************************/

#ifndef __MSA_H
#define __MSA_H

#include "global_variable.h"


extern "C" {
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/utils/basic.h>
    #include <ViennaRNA/io/file_formats_msa.h>

    extern int vrna_file_msa_read(const char * aliFile,char *** Names,char *** seqAlns,char ** Id,char ** structures,unsigned int options);

};


int readMSAsto(string msaFile,vector<string> &Names, vector<string> &seqAlns, string &RNAtype, string &RNAss);
vec_string msafam2prf(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns);
vec_string msa2prf(string querySeq,string msaFile,string ss);
vec_string msafam2prf(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns,string type);

vec_string msafam2prfWithWeight(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns,float weight);
vec_string msafam2prfWithWeight(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns);

#endif
