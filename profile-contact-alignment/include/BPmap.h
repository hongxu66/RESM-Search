/*************************************************************************
    > File Name: align.h
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun May 15 16:19:10 2022
 ************************************************************************/

#ifndef __BPmap_H
#define __BPmap_H

#include "global_variable.h"
#include "buildmap.h"
#include "load_data.h"
#include "msa.h"
#include <omp.h>

static omp_lock_t lock;

//vec_int align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx, mtx_double &p_sco_mtx);

//double Falign(double *sco_mtx, int rows, int cols);
inline double sepw(double sep){if(sep <= 4){return 0.50;}else if(sep == 5){return 0.75;}else{return 1.00;}}

//void mapalignment(Map A, Map B, Options opt,vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO);

//void ini_prf_SCO(mtx_double &P_SCO, double prf_w,Map A, Map B);
//void ini_SCO(double sep_x, double sep_y, mtx_double &SCO,Map &A, Map &B);
//vec_int mod_SCO(double do_it,vec_double gap_a,vec_double gap_b, Map &A, Map &B, double gap_e, mtx_double &SCO, mtx_double &P_SCO);
//void chk (Map &A, Map &B, vec_double gap_a,vec_double gap_b, double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& a2b,mtx_double &P_SCO);
//void map2seqswithmsa(Options &OPT);

inline double exp_fast(double x){
    // WARNING fails if |x| > 1024
    //https://codingforspeed.com/using-faster-exponential-approximation/
    x = 1 + x/1024;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x;
    return x;
}

inline double gaussian(double mean, double stdev, double x){return exp_fast(-pow((x - mean),2)/(2*(pow(stdev,2))));}


//void mapalignment(Map a, Map b, Options opt,vec_double gap_a, vec_double gap_b, mtx_double P_SCO);
//void mapalignment(Map a, Map b, Options opt,vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO);
//void mapalignment(Map a, Map b, Options opt,vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO,string qseq,stringstream &buf);
////double mapalignment(Map a, Map b, Options opt, vec_double gap_a, vec_double gap_b, mtx_double P_SCO, string qseq, string tseq, stringstream &buf);
//void mapalignment(Map a, Map b, Options opt, vec_double gap_a, vec_double gap_b, mtx_double P_SCO, string qseq, stringstream &buf);
//void mapalignment(Map a, Map b, Options opt, vec_double gap_a, vec_double gap_b, mtx_double P_SCO, string qseq, tuple<float,vec_int> &tmpTuple);
//vec_int mapalignment(Map a, Map b, Options opt, vec_double gap_a, vec_double gap_b, string tseq, mtx_double P_SCO);


void CmapSearchDatabase(Options &OPT);
void CmapSearchGenome(Options &OPT);
void CmapAlign(Options &OPT);
void CmapAlignNCBI(Options &OPT);
void CmapAlignNCBI1201(Options &OPT);
void CmapAlignNCBI1202(Options &OPT);
void CmapAli(Options &OPT);
void mulmap2mulseqs(Options &OPT);
void map2seqswithmsa(Options &OPT);
void mulmap2genome(Options &OPT);
void ini_prf_SCO_with_PAM(mtx_double &P_SCO,Map &A, Map &B,int num);
void ini_prf_SCO(mtx_double &P_SCO, double prf_w,Map A, Map B);


class Alignment{
    public:
        Alignment(Map A, Map B, Options OPT,string seqA,string seqB):a(A),b(B),opt(OPT),seqa(seqA),seqb(seqB){};
		Alignment(){};
        ~Alignment(){}
        vec_int mapalignment(vec_double gap_a,vec_double gap_b, mtx_double P_SCO);
        // Return best alignment score; returns 0 when no valid mapping was found.
        double mapalignment(vec_double gap_a, vec_double gap_b, mtx_double P_SCO, stringstream &buf);
//        void mapalignment(vec_double gap_a, vec_double gap_b, mtx_double P_SCO);
        int mapalignment(vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO);
        int mapalignment(vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO,stringstream &buf);
//        void mapalignment(vec_double gap_a, vec_double gap_b, mtx_double P_SCO, stringstream &buf);
//        void mapalignment(vec_double gap_a, vec_double gap_b, mtx_double P_SCO, tuple<float,vec_int> &tmpTuple);

        mtx_int align_label;
        Map a;
        Map b;
        Options opt;
        string nameA;
        string nameB;
        string seqa;
        string seqb;
};

class MSA{
    public:
        std::map<char,int> seqs;
        vector<tuple<string,string,string,string,double>> alignResults;
        vector<vector<string>> clustID;
        void readseqs(vector<string> ids, vector<string> seqs);
        void readCmapResult(string cmapfile);
        void updateAli();

    private:
        string refseq;
        vec_string allseqs;
        vec_string allIDs;

};




int seqsalign(string &seqaid, string &seqa,string &ss,string &seqbid, string &seqb,Options &OPT);
int map2seqsalign(Map &A,string &seqa,string &seqbid, string &seqb,Options &OPT);
vec_int align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx, mtx_double &p_sco_mtx);
double Falign(double *sco_mtx, int rows, int cols);
void ini_SCO(double sep_x, double sep_y, mtx_double &SCO,Map &a, Map &b);
vec_int mod_SCO(double do_it,vec_double gap_a,vec_double gap_b, Map &A, Map &B, double gap_e, mtx_double &SCO, mtx_double &P_SCO);
vec_int mod_SCO(double do_it,vec_double gap_a,vec_double gap_b, Map &A, Map &B, double gap_e, mtx_double &SCO, mtx_double &P_SCO, string seq);
void chk (Map &A, Map &B, vec_double gap_a, vec_double gap_b,double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& a2b,mtx_double &P_SCO);
void chk (Map &A, Map &B, vec_double gap_a, vec_double gap_b,double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& a2b,mtx_double &P_SCO,string qseq);


vec_string a2b2a2m(vec_int a2b, vec_string seq, int modelLen);
vec_string updatealn(vector<vec_int> alla2b,vector<string> seqs,string refseq);
vec_string updatealnnew(vector<vec_int> alla2b,vector<string> seqs,string refseq);
vec_string updatealnnew1(vector<vec_int> alla2b,vector<string> seqs,string refseq,vec_string allIDs);

#endif
