/*************************************************************************
    > File Name: align.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Wed May 25 09:25:39 2022
 ************************************************************************/

#include "align.h"
#include<iostream>
#include<sstream>
using namespace std;

namespace {
clock_t resm_run_start = 0;

int alignment_threads(const Options &opt, int rows, int cols) {
    if (omp_in_parallel()) {
        return 1;
    }
    if ((long long)rows * cols < 40000) {
        return 1;
    }
    return max(1, opt.numThreads);
}
}

void set_resm_run_start(clock_t start_time) {
    resm_run_start = start_time;
}

double get_resm_run_time() {
    if(resm_run_start == 0) {
        return 0;
    }
    return double(clock() - resm_run_start) / CLOCKS_PER_SEC;
}



int Alignment::mapalignment(vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO, stringstream &buf){

//    std::stringstream buf;
    double gap_ext_w = fabs(opt.gap_ext)/fabs(opt.gap_open);

    int max_sep_x = 0;
    int max_sep_y = 0;
    int max_g_e = 0;
    double con_max = -1;
    double gap_max = 0;
    double prf_max = 0;
    vec_int a2b_max;

    int size_a = a.mtx.size();
    int size_b = b.mtx.size();

    // try different sep (sequence seperation difference) penalities
    vec_double sep_x_steps {0,1,2}; // (constant, linear, quadratic)
    vec_double sep_y_steps {1,2,4,8,16,32};
    vec_double gap_e_steps {5,10,100,1000};
    
    // OpenMP parallel scan of parameter space - collapse 3 nested loops
    int threads = alignment_threads(opt, size_a, size_b);
    #pragma omp parallel for collapse(3) schedule(static) num_threads(threads) if(threads > 1)
    for(int sx = 0; sx < sep_x_steps.size(); sx++){
        for(int sy = 0; sy < sep_y_steps.size(); sy++){
            for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){
                double sep_x = sep_x_steps[sx];
                double sep_y = sep_y_steps[sy];
                double gap_e = 1/gap_e_steps[g_e];

                // Get initial score matrix
                mtx_double SCO(size_a,vector<double>(size_b,0));
                ini_SCO(sep_x,sep_y,SCO,a,b);

                // get alignment (a2b mapping) after X iterations
                vec_int a2b = mod_SCO(opt.iter,gap_a,gap_b,a,b,gap_e,SCO,P_SCO,seqa);

                // compute number of contacts/gaps made
                double con_sco = 0;
                double gap_sco = 0;
                double prf_sco = 0;
                chk(a,b,gap_a,gap_b,gap_ext_w,con_sco,gap_sco,prf_sco,a2b,P_SCO,seqa);

                // save if BEST! - use critical section for thread-safe update
                double current_score = con_sco+gap_sco+prf_sco;
                #pragma omp critical(alignment_update)
                {
                    if(current_score > con_max+gap_max+prf_max){
                        max_sep_x = sep_x;
                        max_sep_y = sep_y;
                        max_g_e = g_e;
                        con_max = con_sco;
                        gap_max = gap_sco;
                        prf_max = prf_sco;
                        a2b_max = a2b;
                    }
                }
            }
        }
    }
    // Report the BEST score
    if(a2b_max.empty()){
        buf << "Two maps have not mapped!!!\n";
        return 0;
    }

    double bestScore = con_max + gap_max + prf_max;
    buf << "TMP\t" << max_sep_x << "_" << max_sep_y << "\t" << con_max << "\t" << gap_max << "\t" << prf_max << "\t" << bestScore << "\n";

    int aln_len = 0;for(int ai = 0; ai < size_a; ai++){int bi = a2b_max[ai];if(bi != -1){aln_len++;}}

    vector<vector<int>> maps;
    vector<int> tema,temb;
    if(a2b_max.size()){
        for(int a = 0; a < size_a; a++){int b = a2b_max[a];/*cout << a << " " << b << endl;*/if(b != -1){buf << /*" A\t"*/" " << m2n_a[a] << ":" << m2n_b[b]; tema.push_back(m2n_a[a]); temb.push_back(m2n_b[b]);}}

        if (!tema.empty() && !temb.empty()){
            maps.push_back(tema);
            maps.push_back(temb);

//            alignseqs(opt.seqa,opt.seqb,maps);
        }
        else{
            buf <<"There is no mapping between seqA and seqB!!!" << endl;
        }
    }
    else buf << "Two maps have not mapped!!!\n";
//    return buf;
    return 0;
}

double Alignment::mapalignment(vec_double gap_a, vec_double gap_b, mtx_double P_SCO,  stringstream &buf){

//    std::stringstream buf;
//    buf << tseq << ":test1" << endl;
    double gap_ext_w = fabs(opt.gap_ext)/fabs(opt.gap_open);

    int max_sep_x = 0;
    int max_sep_y = 0;
    int max_g_e = 0;
    double con_max = -1;
    double gap_max = 0;
    double prf_max = 0;
    vec_int a2b_max;

    int size_a = a.mtx.size();
    int size_b = b.mtx.size();

    // try different sep (sequence seperation difference) penalities
    vec_double sep_x_steps {0,1,2}; // (constant, linear, quadratic)
    vec_double sep_y_steps {1,2,4,8,16,32};
    vec_double gap_e_steps {5,10,100,1000};
    
    // OpenMP parallel scan of parameter space - collapse 3 nested loops
    int threads = alignment_threads(opt, size_a, size_b);
    #pragma omp parallel for collapse(3) schedule(static) num_threads(threads) if(threads > 1)
    for(int sx = 0; sx < sep_x_steps.size(); sx++){
        for(int sy = 0; sy < sep_y_steps.size(); sy++){
            for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){
                double sep_x = sep_x_steps[sx];
                double sep_y = sep_y_steps[sy];
                double gap_e = 1/gap_e_steps[g_e];

                // Get initial score matrix
                mtx_double SCO(size_a,vector<double>(size_b,0));
                ini_SCO(sep_x,sep_y,SCO,a,b);

                // get alignment (a2b mapping) after X iterations
                vec_int a2b = mod_SCO(opt.iter,gap_a,gap_b,a,b,gap_e,SCO,P_SCO,seqa);

                // compute number of contacts/gaps made
                double con_sco = 0;
                double gap_sco = 0;
                double prf_sco = 0;
                chk(a,b,gap_a,gap_b,gap_ext_w,con_sco,gap_sco,prf_sco,a2b,P_SCO,seqa);

                // save if BEST! - use critical section for thread-safe update
                double current_score = con_sco+gap_sco+prf_sco;
                #pragma omp critical(alignment_update)
                {
                    if(current_score > con_max+gap_max+prf_max){
                        max_sep_x = sep_x;
                        max_sep_y = sep_y;
                        max_g_e = g_e;
                        con_max = con_sco;
                        gap_max = gap_sco;
                        prf_max = prf_sco;
                        a2b_max = a2b;
                    }
                }
            }
        }
    }
    // Report the BEST score
    if(a2b_max.empty()){
        buf << "Two maps have not mapped!!!\n";
        return 0;
    }

    double bestScore = con_max + gap_max + prf_max;
    buf << "TMP\t" << max_sep_x << "_" << max_sep_y << "\t" << con_max << "\t" << gap_max << "\t" << prf_max << "\t" << bestScore << "\n";

    int aln_len = 0;for(int ai = 0; ai < size_a; ai++){int bi = a2b_max[ai];if(bi != -1){aln_len++;}}

    vector<vector<int>> maps;
    vector<int> tema,temb;
    if(a2b_max.size()){
        for(int i = 0; i < size_a; i++){int j = a2b_max[i];/*cout << a << " " << b << endl;*/if(j != -1){ /*buf << " A\t" " " << a.M2N[i] << ":" << b.M2N[j];*/  tema.push_back(a.M2N[i]); temb.push_back(b.M2N[j]);}}

        if (!tema.empty() && !temb.empty()){
            maps.push_back(tema);
            maps.push_back(temb);
            alignseqs(seqa,seqb,maps);
            buf << ">RESM_Search:" << a.id << ":" << b.id << " SCORE:" << bestScore << " RUN TIME:" << get_resm_run_time() << endl;
//            alignseqs(opt.seqa,opt.seqb,maps);
        }
        else{
            buf <<"There is no mapping between seqA and seqB!!!" << endl;
            return 0;
        }
    }
    else{
        buf << "Two maps have not mapped!!!\n";
        return 0;
    }
    return bestScore;
}

int Alignment::mapalignment(vec_int m2n_a, vec_int m2n_b,vec_double gap_a, vec_double gap_b, mtx_double P_SCO){

    double gap_ext_w = fabs(opt.gap_ext)/fabs(opt.gap_open);

    int max_sep_x = 0;
    int max_sep_y = 0;
    int max_g_e = 0;
    double con_max = -1;
    double gap_max = 0;
    double prf_max = 0;
    vec_int a2b_max;

    int size_a = a.mtx.size();
    int size_b = b.mtx.size();

    // try different sep (sequence seperation difference) penalities
    vec_double sep_x_steps {0,1,2}; // (constant, linear, quadratic)
    vec_double sep_y_steps {1,2,4,8,16,32};
    vec_double gap_e_steps {5,10,100,1000};
    
    // OpenMP parallel scan of parameter space - collapse 3 nested loops
    int threads = alignment_threads(opt, size_a, size_b);
    #pragma omp parallel for collapse(3) schedule(static) num_threads(threads) if(threads > 1)
    for(int sx = 0; sx < sep_x_steps.size(); sx++){
        for(int sy = 0; sy < sep_y_steps.size(); sy++){
            for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){
                double sep_x = sep_x_steps[sx];
                double sep_y = sep_y_steps[sy];
                double gap_e = 1/gap_e_steps[g_e];
            
                // Get initial score matrix
                mtx_double SCO(size_a,vector<double>(size_b,0));
                ini_SCO(sep_x,sep_y,SCO,a,b);         

                // get alignment (a2b mapping) after X iterations
                vec_int a2b = mod_SCO(opt.iter,gap_a,gap_b,a,b,gap_e,SCO,P_SCO);
                
                // compute number of contacts/gaps made
                double con_sco = 0;
                double gap_sco = 0;
                double prf_sco = 0;
                chk(a,b,gap_a,gap_b,gap_ext_w,con_sco,gap_sco,prf_sco,a2b,P_SCO);
              
                // save if BEST! - use critical section for thread-safe update
                double current_score = con_sco+gap_sco+prf_sco;
                #pragma omp critical(alignment_update)
                {
                    if(current_score > con_max+gap_max+prf_max){
                        max_sep_x = sep_x;
                        max_sep_y = sep_y;
                        max_g_e = g_e;
                        con_max = con_sco;
                        gap_max = gap_sco;
                        prf_max = prf_sco;
                        a2b_max = a2b;
                    }
                }
            }
        }
    }
    // Report the BEST score
    if(a2b_max.empty()){
        cout << "Two maps have not mapped!!!\n";
        return 0;
    }
    double bestScore = con_max + gap_max + prf_max;

    int aln_len = 0;for(int ai = 0; ai < size_a; ai++){int bi = a2b_max[ai];if(bi != -1){aln_len++;}}
    
    vector<vector<int>> maps;
    vector<int> tema,temb;
    if(a2b_max.size()){
        for(int a = 0; a < size_a; a++){int b = a2b_max[a];/*cout << a << " " << b << endl;*/if(b != -1){/*cout << " A\t" << m2n_a[a] << ":" << m2n_b[b];*/ tema.push_back(m2n_a[a]); temb.push_back(m2n_b[b]);}}

        if (!tema.empty() && !temb.empty()){
            maps.push_back(tema);
            maps.push_back(temb);
        
            cout << ">RESM_Search:" << a.id << ":" << b.id << " SCORE:" << bestScore << " RUN TIME:" << get_resm_run_time() << endl;
            alignseqs(seqa,seqb,maps);
        }
        else{
            cout <<"There is no mapping between seqA and seqB!!!" << endl;
        }
    }
    else cout << "Two maps have not mapped!!!\n";
    return 0;
}

vec_int Alignment::mapalignment(vec_double gap_a, vec_double gap_b, mtx_double P_SCO){

    double gap_ext_w = fabs(opt.gap_ext)/fabs(opt.gap_open);

    int max_sep_x = 0;
    int max_sep_y = 0;
    int max_g_e = 0;
    double con_max = -1;
    double gap_max = 0;
    double prf_max = 0;
    vec_int a2b_max;

    int size_a = a.mtx.size();
    int size_b = b.mtx.size();

    // try different sep (sequence seperation difference) penalities
    vec_double sep_x_steps {0,1,2}; // (constant, linear, quadratic)
    vec_double sep_y_steps {1,2,4,8,16,32};
    vec_double gap_e_steps {5,10,100,1000};
    
    // OpenMP parallel scan of parameter space - collapse 3 nested loops
    int threads = alignment_threads(opt, size_a, size_b);
    #pragma omp parallel for collapse(3) schedule(static) num_threads(threads) if(threads > 1)
    for(int sx = 0; sx < sep_x_steps.size(); sx++){
        for(int sy = 0; sy < sep_y_steps.size(); sy++){
            for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){
                double sep_x = sep_x_steps[sx];
                double sep_y = sep_y_steps[sy];
                double gap_e = 1/gap_e_steps[g_e];

                // Get initial score matrix
                mtx_double SCO(size_a,vector<double>(size_b,0));
                ini_SCO(sep_x,sep_y,SCO,a,b);

                // get alignment (a2b mapping) after X iterations
                vec_int a2b = mod_SCO(opt.iter,gap_a,gap_b,a,b,gap_e,SCO,P_SCO);

                // compute number of contacts/gaps made
                double con_sco = 0;
                double gap_sco = 0;
                double prf_sco = 0;
                chk(a,b,gap_a,gap_b,gap_ext_w,con_sco,gap_sco,prf_sco,a2b,P_SCO);

                // save if BEST! - use critical section for thread-safe update
                double current_score = con_sco+gap_sco+prf_sco;
                #pragma omp critical(alignment_update)
                {
                    if(current_score > con_max+gap_max+prf_max){
                        max_sep_x = sep_x;
                        max_sep_y = sep_y;
                        max_g_e = g_e;
                        con_max = con_sco;
                        gap_max = gap_sco;
                        prf_max = prf_sco;
                        a2b_max = a2b;
                    }
                }
            }
        }
    }
    // Report the BEST score
    if(a2b_max.empty()){
        cout << "Two maps have not mapped!!!\n";
        vec_int a2b_t;
        return a2b_t;
    }
    double bestScore = con_max + gap_max + prf_max;

    int aln_len = 0;for(int ai = 0; ai < size_a; ai++){int bi = a2b_max[ai];if(bi != -1){aln_len++;}}

    vector<vector<int>> maps;
    vector<int> tema,temb;
    vec_string temseq;
    for(int i=0;i<seqb.length();i++){string tembase(1,seqb[i]);temseq.push_back(tembase);}
//    vec_string temoutseq = a2b2a2m(a2b_max,temseq,size_b);
    if(a2b_max.size()){
        for(int i = 0; i < size_a; i++){int j = a2b_max[i];/*cout << a << " " << b << endl;*/if(j != -1){/*cout << " A\t" << m2n_a[a] << ":" << m2n_b[b];*/ tema.push_back(a.M2N[i]); temb.push_back(b.M2N[j]);}}

        if (!tema.empty() && !temb.empty()){
            maps.push_back(tema);
            maps.push_back(temb);

            cout << ">RESM_Search:" << a.id << ":" << b.id << " SCORE:" << bestScore << " RUN TIME:" << get_resm_run_time() << endl;
            alignseqs(seqa,seqb,maps);
        }
        else{
            cout <<"There is no mapping between seqA and seqB!!!" << endl;
            vec_int a2b_t;
            return a2b_t;
        }

    }
    else {
        cout << "Two maps have not mapped!!!\n";
        vec_int a2b_t;
        return a2b_t;
    }
    return a2b_max;
}


void CmapSearchGenome(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    vec_string bNames, bSeqs;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);
    readfastafile(bNames,bSeqs,OPT.fileb);
	cout << bNames.size() << " : " << bSeqs.size() <<endl;
    for(int j=0;j < aNames.size();j++){
        string seqa;
        string seqaid = aNames[j];
        string seqss = aSS[j];
        seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;


        Map A;
		vec_int m2n_a;
        //vec_int m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.range_a,OPT.prob);
        if(OPT.msa) {
            vec_string msaNames, msaAlns;
            string alnRNAtype,alnRNAss;
            int alnNum = readMSAsto(OPT.msafile,msaNames,msaAlns,alnRNAtype,alnRNAss);

            if(alnNum==0) {cout << "MSAfile is empth!!";}
            m2n_a = A.buildmapswithMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
//            m2n_a = A.buildmapswithMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
        }
        else
            m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.prob);
        A.printout("TMP.map");
        int size_a = A.mtx.size();

        vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}


        for(int i=0;i < bSeqs.size();i++)
        {


            string seqbid = bNames[i];
            string seqb;
			seqb = bSeqs[i];
            transform(seqb.begin(),seqb.end(),seqb.begin(),::toupper);
//            OPT.seqb = seqb;
//            if(count(seqb.begin(),seqb.end(),'X')) continue;
            int seqLength = seqb.length();
            int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;

//                omp_init_lock(&lock);
            omp_set_num_threads(OPT.numThreads);
            #pragma omp parallel for if(K > 1)
            for(int ii = 0;ii < K;ii++){
//                    omp_set_lock(&lock);
                string temseq;
                tuple<float,vec_int> tmpTuple;
                if(ii == K -1 ) temseq = seqb.substr((OPT.lengthE - OPT.overlap)  * ii);
                else
                    temseq = seqb.substr((OPT.lengthE - OPT.overlap) * ii,OPT.lengthE);
                stringstream buf;
                buf << "GenomeWindow:" << i << ":" << ii << ":" << OPT.lengthE << ":" << OPT.overlap << endl;
                Map B;
//                    OPT.seqb = temseq;
                string type = "fasta";
                vec_int m2n_b = B.seq2data(temseq,OPT.sep_cutoff,OPT.prob, type);
                int size_b = B.mtx.size();
                vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

                mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,A,B);}

//                    mapalignment(A,B,OPT,m2n_a,m2n_b,gap_a,gap_b,P_SCO,seqa,buf);
                Alignment Ali(A,B,OPT,seqa,seqb);
                double bestScore = Ali.mapalignment(gap_a,gap_b,P_SCO,buf);
//                    mapalignment(A,B,OPT,gap_a,gap_b,P_SCO,seqa,tmpTuple);
                buf << "Genome:" << ii <<" : " << temseq.length() << "\n";

                #pragma omp critical
//                    allresults.push_back(std::make_tuple(ii,tmpTuple));
                if (bestScore > 0){
                    std::cout << buf.rdbuf() << endl;
                }
            }
        }
        if(OPT.msa) break;
    }
}

void CmapSearchDatabase(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    vec_string bNames, bSeqs;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);
    readfastafile(bNames,bSeqs,OPT.fileb);
	cout << bNames.size() << " : " << bSeqs.size() <<endl;
    for(int j=0;j < aNames.size();j++){
        string seqa;
        string seqaid = aNames[j];
        string seqss = aSS[j];
        seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;
        Map A;
		vec_int m2n_a;
        //vec_int m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.range_a,OPT.prob);
        if(OPT.msa) {
            vec_string msaNames, msaAlns;
            string alnRNAtype,alnRNAss;
            int alnNum = readMSAsto(OPT.msafile,msaNames,msaAlns,alnRNAtype,alnRNAss);

            if(alnNum==0) {cout << "MSAfile is empth!!";}
            m2n_a = A.buildmapswithMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
//            m2n_a = A.buildmapswithMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
        }
        else
            m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.prob);
        A.printout("TMP.map");
        int size_a = A.mtx.size();

        vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}

        omp_set_num_threads(OPT.numThreads);
        #pragma omp parallel for if(bSeqs.size() > 1)
        for(int i=0;i < bSeqs.size();i++)
        {
            string seqbid = bNames[i];
            string seqb;
			seqb = bSeqs[i];
			stringstream buf;
			buf << "Sequence:" << seqbid << ":" << endl;
            transform(seqb.begin(),seqb.end(),seqb.begin(),::toupper);
//            OPT.seqb = seqb;
//            if(count(seqb.begin(),seqb.end(),'X')) continue;
            int seqLength = seqb.length();
            Map B;
            string type = "fasta";
            vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
            int size_b = B.mtx.size();
            vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

            mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,A,B);}

            Alignment Ali(A,B,OPT,seqa,seqb);
            double bestScore = Ali.mapalignment(gap_a,gap_b,P_SCO,buf);

            #pragma omp critical
            if (bestScore > 0 ){
//                cout << ">" << seqbid << "_" << bestScore << endl;
//                cout << seqb << endl;
//                cout << "best-score:" << bestScore << endl;
                std::cout << buf.rdbuf() << endl;
            }
        }
        if(OPT.msa) break;
    }
}


void CmapAlign(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    vec_string bNames, bSeqs;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);
    readfastafile(OPT,bNames,bSeqs,OPT.fileb);
//	cout << bNames.size() << " : " << bSeqs.size() <<endl;


	vector<vec_int> alla2b;
	vec_string allseqs;
	vec_string allIDs;
    for(int j=0;j < aNames.size();j++){
//        string seqa;
        string seqaid = aNames[j];
        string seqss = aSS[j];
        string seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;


        Map A(seqaid);
		vec_int m2n_a;
        //vec_int m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.range_a,OPT.prob);
        if(OPT.msa) {
            vec_string msaNames, msaAlns;
            string alnRNAtype,alnRNAss;
            int alnNum = readMSAsto(OPT.msafile,msaNames,msaAlns,alnRNAtype,alnRNAss);
//			cout << "TEST" << endl;
//			cout << alnNum << endl;
            if(alnNum==0) {cout << "MSAfile is empth!!";}
            m2n_a = A.buildmapswithMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
        }
        else
            m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.prob);
        A.printout("TMP.map");
        int size_a = A.mtx.size();

        vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}

        vec_int tmp;
        for(int t=0;t<seqa.length();t++){
            tmp.push_back(t);
        }
        alla2b.push_back(tmp);
        allseqs.push_back(seqa);
        allIDs.push_back(seqaid);


        omp_set_num_threads(OPT.numThreads);
        #pragma omp parallel for if(bSeqs.size() > 1)
        for(int i=0;i < bSeqs.size();i++)
        {
            string seqbid = bNames[i];
            string seqb;
			seqb = bSeqs[i];
//			cout << seqb.max_size() << endl;
//			cout << seqbid << endl;
//			cout << seqb.length() << endl;
            transform(seqb.begin(),seqb.end(),seqb.begin(),::toupper);
//            OPT.seqb = seqb;
//            if(count(seqb.begin(),seqb.end(),'X')) continue;


//            cout << ">" << seqaid <<":" << seqb.length() << endl;
//            cout << ">" << seqaid << ":" << seqbid << endl;
            int seqLength = seqb.length();
            Map B(seqbid);
            string type = "fasta";
            vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
            int size_b = B.mtx.size();
            vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

            mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
//                mapalignment(A,B,OPT,m2n_a,m2n_b,gap_a,gap_b,P_SCO);
//            mapalignment(A,B,OPT,gap_a,gap_b,seqb, P_SCO);
            Alignment Ali(B,A,OPT,seqb,seqa);
            vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
            #pragma omp critical
            if(a2b_max.size())
            {
            alla2b.push_back(a2b_max);
            allseqs.push_back(seqb);
            allIDs.push_back(seqbid);}
        }
//    cout << alla2b.size() << " " << allseqs.size() << endl;
    vec_string alignedSeqs = updatealnnew1(alla2b,allseqs,seqa,allIDs);

//    vec_string alignedSeqs = updatealnnew(alla2b,allseqs,seqa);
//    cout << "ALIGNEDSEQS" << alignedSeqs.size() <<  endl;
//    for(int i=0;i<alignedSeqs.size();i++){
//        cout << ">" << allIDs[i] << endl;
//        cout << alignedSeqs[i] << endl;
//    }

    }
}




void CmapAlignNCBI(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);

//	cout << bNames.size() << " : " << bSeqs.size() <<endl;


	vector<vec_int> alla2b;
	vec_string allseqs;
	vec_string allIDs;
    for(int j=0;j < aNames.size();j++){
//        string seqa;
        string seqaid = aNames[j];
        string seqss = aSS[j];
        string seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;


        Map A(seqaid);
		vec_int m2n_a;
        //vec_int m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.range_a,OPT.prob);
        if(OPT.msa) {
            vec_string msaNames, msaAlns;
            string alnRNAtype,alnRNAss;
            int alnNum = readMSAsto(OPT.msafile,msaNames,msaAlns,alnRNAtype,alnRNAss);
//			cout << "TEST" << endl;
//			cout << alnNum << endl;
            if(alnNum==0) {cout << "MSAfile is empth!!";}
            m2n_a = A.buildmapswithMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
        }
        else
            m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.prob);
        A.printout("TMP.map");
        int size_a = A.mtx.size();

        vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}

        vec_int tmp;
        for(int t=0;t<seqa.length();t++){
            tmp.push_back(t);
        }
        alla2b.push_back(tmp);
        allseqs.push_back(seqa);
        allIDs.push_back(seqaid);

        vec_string bNames, bSeqs;
        ifstream in(OPT.fileb.c_str());
        if(!in)
        {
            std::cout << "fastafile is empty!!!" << std::endl;
            exit(1);
        }
        else{
            std::string line, name, seq;
            while( std::getline( in, line ) ){
                if (line.empty()) continue;
                if( line[0] == '>' ){
                    if( !name.empty() ){
                        if(seq.length() >= OPT.lengthC){
                        int seqLength = seq.length();
                        int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;
                        for(int ii=0;ii<K;ii++){
                            string temseq;
                            string newname = name + ":" + to_string(ii);
                            if(ii == K -1 ) temseq = seq.substr((OPT.lengthE - OPT.overlap)  * ii);
                            else
                                temseq = seq.substr((OPT.lengthE - OPT.overlap) * ii,OPT.lengthE);
                            bNames.push_back(newname);
                            bSeqs.push_back(replace_all(temseq,"T","U"));
                        }
                    }
                    else{
                        bNames.push_back(name);
                        bSeqs.push_back(replace_all(seq,"T","U"));
//                        name.clear();
                    }
                        name.clear();
                        if(bSeqs.size() >= 320000){
                            omp_set_num_threads(OPT.numThreads);
                            #pragma omp parallel for if(bSeqs.size() > 1)
                            for(int i=0;i < bSeqs.size();i++)
                            {
                                string seqbid = bNames[i];
                                string seqb;
                                seqb = bSeqs[i];
                                transform(seqb.begin(),seqb.end(),seqb.begin(),::toupper);
                                int seqLength = seqb.length();
                                Map B(seqbid);
                                string type = "fasta";
                                vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
                                int size_b = B.mtx.size();
                                vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

                                mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
                                Alignment Ali(B,A,OPT,seqb,seqa);
                                vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
                                #pragma omp critical
                                cout  << endl;
                            }
                            vector <string>().swap(bNames);
                            vector <string>().swap(bSeqs);
                        }
                    }
                    if( !line.empty() ){
                        name = line.substr(1);
                    }
                    seq.clear();
                }
                else if( !name.empty() ){
                    if( line.find(' ') != std::string::npos ){
                        name.clear();
                        seq.clear();
                    } else {
                        seq += line;
                    }
                }
            }
            if( !name.empty() ){
                //std::cout << name << " : " << seq << std::endl;
                if(seq.length() >= OPT.lengthC){
                        int seqLength = seq.length();
                        int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;
                        for(int ii=0;ii<K;ii++){
                            string temseq;
                            string newname = name + ":" + to_string(ii);
                            if(ii == K -1 ) temseq = seq.substr((OPT.lengthE - OPT.overlap)  * ii);
                            else
                                temseq = seq.substr((OPT.lengthE - OPT.overlap) * ii,OPT.lengthE);
                            bNames.push_back(newname);
                            bSeqs.push_back(replace_all(temseq,"T","U"));
                        }
                    }
                    else{
                        bNames.push_back(name);
                        bSeqs.push_back(replace_all(seq,"T","U"));
//                        name.clear();
                    }
    //			cout << seq << endl;
            }
        }


        omp_set_num_threads(OPT.numThreads);
        #pragma omp parallel for if(bSeqs.size() > 1)
        for(int i=0;i < bSeqs.size();i++)
        {
            string seqbid = bNames[i];
            string seqb;
            seqb = bSeqs[i];
            transform(seqb.begin(),seqb.end(),seqb.begin(),::toupper);
            int seqLength = seqb.length();
            Map B(seqbid);
            string type = "fasta";
            vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
            int size_b = B.mtx.size();
            vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

            mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
    //                mapalignment(A,B,OPT,m2n_a,m2n_b,gap_a,gap_b,P_SCO);
    //            mapalignment(A,B,OPT,gap_a,gap_b,seqb, P_SCO);
            Alignment Ali(B,A,OPT,seqb,seqa);
            vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
            #pragma omp critical
            cout  << endl;
//            if(a2b_max.size())
//            {
//            alla2b.push_back(a2b_max);
//            allseqs.push_back(seqb);
//            allIDs.push_back(seqbid);}
        }
//    cout << alla2b.size() << " " << allseqs.size() << endl;
//    vec_string alignedSeqs = updatealnnew1(alla2b,allseqs,seqa,allIDs);

//    vec_string alignedSeqs = updatealnnew(alla2b,allseqs,seqa);
//    cout << "ALIGNEDSEQS" << alignedSeqs.size() <<  endl;
//    for(int i=0;i<alignedSeqs.size();i++){
//        cout << ">" << allIDs[i] << endl;
//        cout << alignedSeqs[i] << endl;
//    }

    }
}

int map2seqsalign(Map &A,string &seqa, string &seqbid, string &seqb,Options &OPT){
    string type = "fasta";
//    Map A(seqaid);
//    vec_int m2n_a = A.buildmaps(seqa,ss,OPT.sep_cutoff,OPT.prob);
    int size_a = A.mtx.size();
    vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}

    int seqLength = seqb.length();
    if(seqLength > OPT.lengthC) {
        int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;
        for(int i=0;i<K;i++){
            string temseq;
            string temseqid;
            if(i == K -1 )
                temseq = seqb.substr((OPT.lengthE - OPT.overlap)  * i);
            else
                temseq = seqb.substr((OPT.lengthE - OPT.overlap) * i,OPT.lengthE);
                temseqid = seqbid + ":" + std::to_string(i);
                Map B(temseqid);
                vec_int m2n_b = B.seq2data(temseq,OPT.sep_cutoff,OPT.prob, type);
                int size_b = B.mtx.size();
                vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}
                mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
                Alignment Ali(B,A,OPT,temseq,seqa);
                vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
        }
        return 0;
    }

    Map B(seqbid);

    vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
    int size_b = B.mtx.size();
    vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}
    mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
    Alignment Ali(B,A,OPT,seqb,seqa);
    vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
    return 0;
}

int seqsalign(string &seqaid, string &seqa,string &ss,string &seqbid, string &seqb,Options &OPT){
    string type = "fasta";
    Map A(seqaid);
    vec_int m2n_a = A.buildmaps(seqa,ss,OPT.sep_cutoff,OPT.prob);
    int size_a = A.mtx.size();
    vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}

    int seqLength = seqb.length();
    if(seqLength > OPT.lengthC) {
        int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;
        for(int i=0;i<K;i++){
            string temseq;
            string temseqid;
            if(i == K -1 )
                temseq = seqb.substr((OPT.lengthE - OPT.overlap)  * i);
            else
                temseq = seqb.substr((OPT.lengthE - OPT.overlap) * i,OPT.lengthE);
                temseqid = seqbid + ":" + std::to_string(i);
                Map B(temseqid);
                vec_int m2n_b = B.seq2data(temseq,OPT.sep_cutoff,OPT.prob, type);
                int size_b = B.mtx.size();
                vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}
                mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
                Alignment Ali(B,A,OPT,temseq,seqa);
                vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
        }
        return 0;
    }

    Map B(seqbid);

    vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
    int size_b = B.mtx.size();
    vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}
    mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,B,A);}
    Alignment Ali(B,A,OPT,seqb,seqa);
    vec_int a2b_max = Ali.mapalignment(gap_b,gap_a, P_SCO);
    return 0;
}


void CmapAlignNCBI1202(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);

    vec_string bNames, bSeqs;
    readfastafile(bNames,bSeqs,OPT.fileb);

    for(int j=0;j < aNames.size();j++){
        string seqaid = aNames[j];
        string seqss = aSS[j];
        string seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;
        string type = "fasta";
        Map A(seqaid);
        vec_int m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.prob);
//    int size_a = A.mtx.size();
//    vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}
        omp_set_num_threads(OPT.numThreads);
        #pragma omp parallel for if(bSeqs.size() > 1)
        for(int i=0;i<bSeqs.size();i++){
            string seqbid = bNames[i];
            string seqb = bSeqs[i];
            if (count(seqb.begin(),seqb.end(),'X')) continue;
            int tmp = map2seqsalign(A,seqa,seqbid,seqb,OPT);
            #pragma omp critical
            cout << endl;
        }
    }
}
void CmapAlignNCBI1201(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);

    vec_string bNames, bSeqs;
    readfastafile(bNames,bSeqs,OPT.fileb);

    for(int j=0;j < aNames.size();j++){
        string seqaid = aNames[j];
        string seqss = aSS[j];
        string seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;
        omp_set_num_threads(OPT.numThreads);
        #pragma omp parallel for if(bSeqs.size() > 1)
        for(int i=0;i<bSeqs.size();i++){
            string seqbid = bNames[i];
            string seqb = bSeqs[i];
            if (count(seqb.begin(),seqb.end(),'X')) continue;
            int tmp = seqsalign(seqaid,seqa,seqss,seqbid,seqb,OPT);
            #pragma omp critical
            cout << endl;
        }
    }
}

void CmapAli(Options &OPT){
    vec_string aNames, aSeqs, aSS;
    vec_string bNames, bSeqs;
    readdbnfile(aNames,aSeqs,aSS,OPT.filea);

    vec_string bmsaNames,bmsaAlns;
    readfastafile(bmsaNames,bmsaAlns,OPT.msafileb);

    readfastafile(bNames,bSeqs,OPT.fileb);
	cout << bNames.size() << " : " << bSeqs.size() <<endl;
    for(int j=0;j < aNames.size();j++){
        string seqa;
        string seqaid = aNames[j];
        string seqss = aSS[j];
        seqa = aSeqs[j];
        transform(seqa.begin(),seqa.end(),seqa.begin(),::toupper);
//        OPT.seqa = seqa;
        if (count(seqa.begin(),seqa.end(),'X')) continue;


        Map A;
		vec_int m2n_a;
        //vec_int m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.range_a,OPT.prob);
        if(OPT.msa) {
            vec_string msaNames, msaAlns;
            string alnRNAtype,alnRNAss;
//            int alnNum = readMSAsto(OPT.msafile,msaNames,msaAlns,alnRNAtype,alnRNAss);
            readfastafile(msaNames,msaAlns,OPT.msafile);
            int alnNum = msaAlns.size();
			cout << "TEST" << endl;
			cout << alnNum << endl;
            if(alnNum==0) {cout << "MSAfile is empth!!";}
            m2n_a = A.buildmapsMSA(seqa,seqaid,seqss,msaNames,msaAlns,OPT.sep_cutoff);
        }
        else
            m2n_a = A.buildmaps(seqa,seqss,OPT.sep_cutoff,OPT.prob);
        A.printout("TMP.map");
        int size_a = A.mtx.size();

        vec_double gap_a(size_a,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_a,A.ss,OPT.gap_ss_w);}


        for(int i=0;i < bSeqs.size();i++)
        {


            string seqbid = bNames[i];
            string seqb;
			seqb = bSeqs[i];
//			cout << seqb.max_size() << endl;
//			cout << seqbid << endl;
//			cout << seqb.length() << endl;
            transform(seqb.begin(),seqb.end(),seqb.begin(),::toupper);
//            OPT.seqb = seqb;
//            if(count(seqb.begin(),seqb.end(),'X')) continue;


            cout << seqaid <<":" << seqb.length() << endl;
            int seqLength = seqb.length();
            if(seqLength <= OPT.lengthC) {
                Map B;
                string type = "fasta";
//                vec_int m2n_b = B.seq2data(seqb,OPT.sep_cutoff,OPT.prob, type);
                vec_int m2n_b = B.buildmapsMSA(seqb,seqbid,bmsaNames,bmsaAlns,OPT.sep_cutoff);
                int size_b = B.mtx.size();
                vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

                mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,A,B);}
//                mapalignment(A,B,OPT,m2n_a,m2n_b,gap_a,gap_b,P_SCO);
                Alignment Ali(A,B,OPT,seqa,seqb);
                Ali.mapalignment(gap_a,gap_b,P_SCO);
            }
            else{
                int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;
                string temseq;
                omp_set_num_threads(OPT.numThreads);
                #pragma omp parallel for if(K > 1)
                for(int ii = 0;ii < K;ii++){
                    if(ii == K -1 ) temseq = seqb.substr((OPT.lengthE - OPT.overlap)  * ii);
                    else
                        temseq = seqb.substr((OPT.lengthE - OPT.overlap) * ii,OPT.lengthE);
                    Map B;
                    string type = "fasta";
                    vec_int m2n_b = B.seq2data(temseq,OPT.sep_cutoff,OPT.prob, type);
                    int size_b = B.mtx.size();
                    vec_double gap_b(size_b,OPT.gap_open);if(OPT.use_gap_ss == true){mod_gap(gap_b,B.ss,OPT.gap_ss_w);}

                    mtx_double P_SCO;if(OPT.use_prf == true){ini_prf_SCO(P_SCO,OPT.prf_w,A,B);}
                    stringstream buf;
//                    mapalignment(A,B,OPT,m2n_a,m2n_b,gap_a,gap_b,P_SCO,seqa,buf);
                    Alignment Ali(A,B,OPT,seqa,temseq);
                    Ali.mapalignment(gap_a,gap_b,P_SCO,buf);
                    buf << "\nGenome:" << ii << "\n";
                    #pragma omp critical
                    std::cout << buf.rdbuf();
                }
            }
        }
        if(OPT.msa) break;
    }
}


vec_int align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx, mtx_double &p_sco_mtx,string qseq){
    // LOCAL_ALIGN
    // Start    0
    // [A]lign  1
    // [D]own   2
    // [R]ight  3
    
    double max_sco = 0;
    int rows = sco_mtx.size();
    if(rows == 0 || sco_mtx[0].empty()){
        return vec_int();
    }
    int cols = sco_mtx[0].size();
    
    bool add_prf = false;if(p_sco_mtx.size() == (size_t)rows && (rows == 0 || p_sco_mtx[0].size() == (size_t)cols)){add_prf = true;}
    
    vec_int a2b(rows,-1);
    
    int stride = cols + 1;
    size_t cells = (size_t)(rows + 1) * stride;
    static thread_local vector<double> sco;
    static thread_local vector<unsigned char> label;
    if(sco.size() < cells){sco.resize(cells);}
    if(label.size() < cells){label.resize(cells);}
    fill(sco.begin(), sco.begin() + cells, 0.0);
    fill(label.begin(), label.begin() + cells, 0);
    
    int max_i = 0;int max_j = 0;
    for (int i = 1; i <= rows; i++){
        
        for (int j = 1; j <= cols; j++){
            size_t ij = (size_t)i * stride + j;
            double A = sco[(size_t)(i-1) * stride + (j-1)] + sco_mtx[i-1][j-1]; if(add_prf == true){A += p_sco_mtx[i-1][j-1];}
            double D = sco[(size_t)(i-1) * stride + j];
            double R = sco[(size_t)i * stride + (j-1)];
            
            if(label[(size_t)(i-1) * stride + j] == 1){D += gap_b[j-1];}else{D += gap_b[j-1] * gap_e;}
            if(label[(size_t)i * stride + (j-1)] == 1){R += gap_a[i-1];}else{R += gap_a[i-1] * gap_e;}
            
            if(A <= 0 and D <= 0 and R <= 0){label[ij] = 0;sco[ij] = 0;}
            else{
                if(A >= R){if(A >= D){label[ij] = 1;sco[ij] = A;}else{label[ij] = 2;sco[ij] = D;}}
                else{if(R >= D){label[ij] = 3;sco[ij] = R;}else{label[ij] = 2;sco[ij] = D;}}
                if(sco[ij] > max_sco){max_i = i;max_j = j;max_sco = sco[ij];}
            }
        }
    }
    int i = max_i;int j = max_j;
    while(1){
        unsigned char move = label[(size_t)i * stride + j];
        if(move == 0){break;}
        else if(move == 1){a2b[i-1] = j-1;i--;j--;}
        else if(move == 2){i--;}
        else if(move == 3){j--;}
    }
    return(a2b);
}

vec_int align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx, mtx_double &p_sco_mtx){
    // LOCAL_ALIGN
    // Start    0
    // [A]lign  1
    // [D]own   2
    // [R]ight  3
    
    double max_sco = 0;
    int rows = sco_mtx.size();
    if(rows == 0 || sco_mtx[0].empty()){
        return vec_int();
    }
    int cols = sco_mtx[0].size();
    
    bool add_prf = false;if(p_sco_mtx.size() == (size_t)rows && (rows == 0 || p_sco_mtx[0].size() == (size_t)cols)){add_prf = true;}
    
    vec_int a2b(rows,-1);
    
    int stride = cols + 1;
    size_t cells = (size_t)(rows + 1) * stride;
    static thread_local vector<double> sco;
    static thread_local vector<unsigned char> label;
    if(sco.size() < cells){sco.resize(cells);}
    if(label.size() < cells){label.resize(cells);}
    fill(sco.begin(), sco.begin() + cells, 0.0);
    fill(label.begin(), label.begin() + cells, 0);
    
    int max_i = 0;int max_j = 0;
    for (int i = 1; i <= rows; i++){
        
        for (int j = 1; j <= cols; j++){
            size_t ij = (size_t)i * stride + j;
            double A = sco[(size_t)(i-1) * stride + (j-1)] + sco_mtx[i-1][j-1]; if(add_prf == true){A += p_sco_mtx[i-1][j-1];}
            double D = sco[(size_t)(i-1) * stride + j];
            double R = sco[(size_t)i * stride + (j-1)];
            
            if(label[(size_t)(i-1) * stride + j] == 1){D += gap_b[j-1];}else{D += gap_b[j-1] * gap_e;}
            if(label[(size_t)i * stride + (j-1)] == 1){R += gap_a[i-1];}else{R += gap_a[i-1] * gap_e;}
            
            if(A <= 0 and D <= 0 and R <= 0){label[ij] = 0;sco[ij] = 0;}
            else{
                if(A >= R){if(A >= D){label[ij] = 1;sco[ij] = A;}else{label[ij] = 2;sco[ij] = D;}}
                else{if(R >= D){label[ij] = 3;sco[ij] = R;}else{label[ij] = 2;sco[ij] = D;}}
                if(sco[ij] > max_sco){max_i = i;max_j = j;max_sco = sco[ij];}
            }
        }
    }
    int i = max_i;int j = max_j;
    while(1){
        unsigned char move = label[(size_t)i * stride + j];
        if(move == 0){break;}
        else if(move == 1){a2b[i-1] = j-1;i--;j--;}
        else if(move == 2){i--;}
        else if(move == 3){j--;}
    }
    return(a2b);
}
double Falign(double *sco_mtx, int rows, int cols){
    if(rows <= 0 || cols <= 0){
        return 0;
    }
    double max_sco = 0;
    static thread_local vector<double> sco;
    if(sco.size() < (size_t)cols + 1){sco.resize((size_t)cols + 1);}
    fill(sco.begin(), sco.begin() + cols + 1, 0.0);
    for (int i = 1; i <= rows; i++){
        double diag = 0;
        for (int j = 1; j <= cols; j++){
            double up = sco[j];
            double A = diag + sco_mtx[(i-1)*cols+(j-1)];
            double D = up;
            double R = sco[j-1];
            diag = up;
            
            if(A >= R){if(A >= D){sco[j] = A;}else{sco[j] = D;}}
            else{if(R >= D){sco[j] = R;}else{sco[j] = D;}}
            
            if(sco[j] > max_sco){max_sco = sco[j];}
        }
    }
    return(max_sco);
}


void ini_SCO(double sep_x, double sep_y, mtx_double &SCO,Map &a, Map &b){
    // Get initial score matrix
    // Pre-allocate buffer to avoid repeated VLA allocation - thread local storage
    static thread_local std::vector<double> M_buffer(65536);  // 256x256 max
    
    for(int i=0; i < a.vec.size(); i++){ // go through columns (vec_a) in map_a that has contacts
        int ai = a.vec[i];
        for(int j=0; j < b.vec.size(); j++){ // go through columns (vec_b) in map_b that has contacts
            int bi = b.vec[j];
            int A[2] = {(int)a.vec_div[ai],(int)(a.vec_i[ai].size()-a.vec_div[ai])};
            int B[2] = {(int)b.vec_div[bi],(int)(b.vec_i[bi].size()-b.vec_div[bi])};
            for(int k=0; k <= 1; k++){ // left and right of diagonal
                if(A[k] > 0 and B[k] > 0){
                    int M_size = A[k]*B[k];
                    // Resize buffer if needed
                    if(M_size > M_buffer.size()){
                        M_buffer.resize(M_size * 2);
                    }
                    
                    for(int n=0; n < A[k]; n++){
                        int nn = n; if(k == 1){nn += a.vec_div[ai];}
                        int aj = a.vec_i[ai][nn];
                        int sep_a = abs(ai-aj);
                        for(int m=0; m < B[k]; m++){
                            int mm = m; if(k == 1){mm += b.vec_div[bi];}
                            int bj = b.vec_i[bi][mm];
                            int sep_b = abs(bi-bj);
                            int sep_D = abs(sep_a-sep_b);
                            double sep_M = min(sep_a,sep_b);
                            //double sep_std = sep_y*(1+pow(sep_M-2,sep_x));
                            double sep_std = sep_y*(1+pow(sep_M-1,sep_x));
                            if(sep_D/sep_std < 6){
                                M_buffer[n*B[k]+m] = a.mtx[ai][aj] * b.mtx[bi][bj] * sepw(sep_M) * gaussian(0,sep_std,sep_D);
                            }else{M_buffer[n*B[k]+m] = 0;}
                        }
                    }

                    SCO[ai][bi] += Falign(M_buffer.data(),A[k],B[k]);
                }
            }
        }
    }
}

vec_int mod_SCO(double do_it,vec_double gap_a,vec_double gap_b, Map &A, Map &B, double gap_e, mtx_double &SCO, mtx_double &P_SCO,string qseq){
    // iterate
    vec_int a2b_tmp;
    for(int it=0; it < do_it; it++)
    {
        // align
//        a2b_tmp = align(gap_a,gap_b,gap_e,SCO,P_SCO);
        a2b_tmp = align(gap_a,gap_b,gap_e,SCO,P_SCO,qseq);
        
        // update similarity matrix
        double IT = (double)it + 1;
        double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
        for(int a=0; a < A.vec.size(); a++){ // go through columns (vec_a) in map_a that has contacts
            int ai = A.vec[a];
            for(int b=0; b < B.vec.size(); b++){ // go through columns (vec_b) in map_b that has contacts
                int bi = B.vec[b];
                double sco_contact = 0;
                for(int n=0; n < A.vec_i[ai].size(); n++){ // go through contacts in vec_a
                    int aj = A.vec_i[ai][n];
                    int bj = a2b_tmp[aj]; // get mapping
                    if(bj != -1){ // if mapping exists
                        if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
                            double sep_M = min(abs(ai-aj),abs(bi-bj));
                            sco_contact += A.mtx[ai][aj] * B.mtx[bi][bj] * sepw(sep_M);
                        }
                    }
                }
                SCO[ai][bi] = s1*SCO[ai][bi] + s2*sco_contact;
            }
        }
    }
    return(a2b_tmp);
}


vec_int mod_SCO(double do_it,vec_double gap_a,vec_double gap_b, Map &A, Map &B, double gap_e, mtx_double &SCO, mtx_double &P_SCO){
    // iterate
    vec_int a2b_tmp;
    for(int it=0; it < do_it; it++)
    {
        // align
        a2b_tmp = align(gap_a,gap_b,gap_e,SCO,P_SCO);
        
        // update similarity matrix
        double IT = (double)it + 1;
        double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
        for(int a=0; a < A.vec.size(); a++){ // go through columns (vec_a) in map_a that has contacts
            int ai = A.vec[a];
            for(int b=0; b < B.vec.size(); b++){ // go through columns (vec_b) in map_b that has contacts
                int bi = B.vec[b];
                double sco_contact = 0;
                for(int n=0; n < A.vec_i[ai].size(); n++){ // go through contacts in vec_a
                    int aj = A.vec_i[ai][n];
                    int bj = a2b_tmp[aj]; // get mapping
                    if(bj != -1){ // if mapping exists
                        if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
                            double sep_M = min(abs(ai-aj),abs(bi-bj));
                            sco_contact += A.mtx[ai][aj] * B.mtx[bi][bj] * sepw(sep_M);
                        }
                    }
                }
                SCO[ai][bi] = s1*SCO[ai][bi] + s2*sco_contact;
            }
        }
    }
    return(a2b_tmp);
}

void chk (Map &A, Map &B, vec_double gap_a, vec_double gap_b,double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& a2b,mtx_double &P_SCO,string qseq){
    
    int size_a = A.mtx.size();
    if(a2b.size() < (size_t)size_a){
        return;
    }
    bool use_prf = false; if(P_SCO.size() == (size_t)size_a && (size_a == 0 || P_SCO[0].size() == B.mtx.size())){use_prf = true;}
    
    int a = 0;int b = 0;
    for(int ai = 0; ai < size_a; ai++){
        int bi = a2b[ai];
        if(bi != -1){
            if(use_prf == true){prf_sco += P_SCO[ai][bi];}
            if(a > 0){ // compute number of gaps
                double num_gap_a = ((ai-a)-1); if(qseq[a] == 'N') ; else {if(num_gap_a > 0){gap_sco += gap_a[ai] + gap_a[ai] * gap_e_w * (num_gap_a-1);}}
                double num_gap_b = ((bi-b)-1); if(num_gap_b > 0){gap_sco += gap_b[bi] + gap_b[bi] * gap_e_w * (num_gap_b-1);}
            }
            for(int m=0; m < A.vec_div[ai]; m++){ // compute number of contacts
                int aj = A.vec_i[ai][m];
                int bj = a2b[aj];
                if(bj != -1){
                    double sep_M = min(abs(ai-aj),abs(bi-bj));
                    con_sco += A.mtx[ai][aj] * B.mtx[bi][bj] * sepw(sep_M);
                }
            }
            a = ai;b = bi;
        }
    }
    gap_sco /= 2;
}

void chk (Map &A, Map &B, vec_double gap_a, vec_double gap_b,double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& a2b,mtx_double &P_SCO){
    
    int size_a = A.mtx.size();
    if(a2b.size() < (size_t)size_a){
        return;
    }
    bool use_prf = false; if(P_SCO.size() == (size_t)size_a && (size_a == 0 || P_SCO[0].size() == B.mtx.size())){use_prf = true;}
    
    int a = 0;int b = 0;
    for(int ai = 0; ai < size_a; ai++){
        int bi = a2b[ai];
        if(bi != -1){
            if(use_prf == true){prf_sco += P_SCO[ai][bi];}
            if(a > 0){ // compute number of gaps
                double num_gap_a = ((ai-a)-1); if(num_gap_a > 0){gap_sco += gap_a[ai] + gap_a[ai] * gap_e_w * (num_gap_a-1);}
                double num_gap_b = ((bi-b)-1); if(num_gap_b > 0){gap_sco += gap_b[bi] + gap_b[bi] * gap_e_w * (num_gap_b-1);}
            }
            for(int m=0; m < A.vec_div[ai]; m++){ // compute number of contacts
                int aj = A.vec_i[ai][m];
                int bj = a2b[aj];
                if(bj != -1){
                    double sep_M = min(abs(ai-aj),abs(bi-bj));
                    con_sco += A.mtx[ai][aj] * B.mtx[bi][bj] * sepw(sep_M);
                }
            }
            a = ai;b = bi;
        }
    }
    gap_sco /= 2;
}

// compute profile  matrix

void ini_prf_SCO_with_PAM(mtx_double &P_SCO,Map &A, Map &B,int num){
    int size_a = A.prf.size();
    int size_b = B.prf.size();
    
    P_SCO.resize(size_a,vector<double>(size_b,0));

    Matrix4d PAMatrix = PAM(num);
    for(int i=0;i<size_a;i++)
        for(int j=0;j<size_b;j++){
            if(A.aa[i] == 'N' or B.aa[j] == 'N') P_SCO[i][j] = 0;
            else{
                int k1,k2;
                if (A.aa[i] == 'A') k1 = 0;
                else if(A.aa[i] == 'U') k1 = 1;
                else if(A.aa[i] == 'C') k1 = 2;
                else if(A.aa[i] == 'G') k1 = 3;
                else {cout << "BASE ERROR" << A.aa[i] << endl; exit(1);}
                if (B.aa[j] == 'A') k2 = 0;
                else if(B.aa[j] == 'U') k2 = 1;
                else if(B.aa[j] == 'C') k2 = 2;
                else if(B.aa[j] == 'G') k2 = 3;
                else {cout << "BASE ERROR" << B.aa[j] << endl; exit(1);}
                P_SCO[i][j] = PAMatrix(k1,k2);
            }
        }
}


void ini_prf_SCO(mtx_double &P_SCO, double prf_w,Map A, Map B){
    int size_a = A.prf.size();
    int size_b = B.prf.size();
    
    //cout << "size_a:size_b " << size_a << ":" << size_b << endl;
    P_SCO.resize(size_a,vector<double>(size_b,0));
    if(size_a == 0 || size_b == 0 || A.prf.empty() || B.prf.empty() || A.prf[0].empty()){
        return;
    }
    
    // compute background frequencies
    //vec_double pb(20,0); int prf_size = prf_a[0].size(); // modified by HX at 2022.02.18
    vec_double pb(4,0); int prf_size = A.prf[0].size(); // modified by HX at 2022.02.18
    double pb_size = 0;
    for(int ai = 0; ai < size_a; ai++)
    {
        if(A.aa[ai] != 'N'){ // ignore positions that have no identity
            for(int p=0; p < prf_size; p++){pb[p] += A.prf[ai][p];}
            pb_size += 1;
        }
    }
    for(int bi = 0; bi < size_b; bi++)
    {
        if(B.aa[bi] != 'N'){ // ignore positions that have no identity
            for(int p=0; p < prf_size; p++){pb[p] += B.prf[bi][p];}
            pb_size += 1;
        }
    }
    if(pb_size == 0){
        return;
    }
    for (int i=0; i < size_a; i++){
        for (int j=0; j < size_b; j++){
            if(A.aa[i] == 'N' or B.aa[j] == 'N'){P_SCO[i][j] = 0;} // if no identity, return score of 0
            else{
                // profile comparison calculation, similar to HHsuite from Soeding.
                double tmp_sco = 0;
                for(int p=0; p < prf_size; p++){
                    if(pb[p] > 0){
                        tmp_sco += (A.prf[i][p]*B.prf[j][p])/(pb[p]/pb_size);
                    }
                }
                //P_SCO[i][j] = log2(tmp_sco)/5 * prf_w;
                P_SCO[i][j] = tmp_sco > 0 ? log2(tmp_sco)/4 * prf_w : 0;
            }
        }
    } 
}

vec_string a2b2a2m(vec_int a2b, vec_string seq, int modelLen){
//    cout << "Start a2b2a2m" << endl;
//    for(int i =0;i<a2b.size();i++)
//        cout << a2b[i] << " ";
//    cout << endl;
//    for(int i=0;i<seq.size();i++)
//        cout << seq[i] ;
//     cout << endl;
    int last_idx = -1;
    vec_int a2b0b = a2b;
    int inslen = 0;

    for(int i=0;i<a2b0b.size();i++){
        if (a2b0b[i] == -1){
            for(int j=i+1;j<a2b0b.size();j++){
                if(a2b0b[j]>last_idx)
                    a2b0b[j] += 1;
            }
            a2b0b[i] = last_idx + 1;
            inslen++;
        }
        last_idx = a2b0b[i];
    }
//    cout << "TEST1" << endl;
    vec_string seqa2m(modelLen+inslen,"-");
//    for(int i=0;i<modelLen+inslen;i++){
//        seqa2m.push_back("-");
//    }
//    cout << "TEST1" << endl;
    for(int i=0;i<a2b0b.size();i++){
//        cout << a2b0b[i] << " ";
        string tem = seq[i];
        if(a2b[i] == -1){
            transform(tem.begin(),tem.end(),tem.begin(),::tolower);
            seqa2m[a2b0b[i]] = tem;
        }
        else{
            transform(tem.begin(),tem.end(),tem.begin(),::toupper);
            seqa2m[a2b0b[i]] = tem;
        }
    }
//    cout << "TEST1" << endl;
//    cout << "Alignment:" ;
//    for(int i=0;i<seqa2m.size();i++)
//        cout << seqa2m[i];
//    cout << endl;
    return seqa2m;
}


vec_string updatealnnew(vector<vec_int> alla2b,vector<string> seqs,string refseq){
    vec_string alignedSeqs;
    vector<vec_int> alla2b0b = alla2b;
    int maxLen = 0;
    int seqLen = 0;
    for(int m=0;m<alla2b.size();m++){
        if (alla2b[m].size() > seqLen) seqLen = alla2b[m].size();
        int inslen = 0;
        for(int n=0;n<alla2b[m].size();n++){
            if(alla2b[m][n] == -1) inslen++;
        }
        if(inslen > maxLen) maxLen = inslen;
    }

    for(int m=0;m<alla2b.size();m++){
        for(int n=0;n<alla2b[m].size();n++){
            cout << alla2b[m][n] << " ";
        }
        cout << endl;
    }

    cout << maxLen << endl;
    vec_int tmpList;
    int last_idx = -1;
    int flag = 0;
    for(int i=0;i<seqLen;i++){
        flag = 0;
        for(int j=0;j<alla2b.size();j++){
            cout << alla2b[j][i] << " ";
            if(i>=alla2b[j].size()) continue;
            if(alla2b[j][i] == -1) {flag=1;break;}
        }
        cout << " TEST " << flag;
        if(flag) {
            for(int m=0;m<alla2b.size();m++){
                if (i>=alla2b[m].size()) continue;
                for(int n=i+1;n<alla2b[m].size();n++){
                    if(alla2b[m][n] > last_idx) alla2b[m][n] += 1;
                }
            }
            for(int k=0;k<alla2b.size();k++){
                if (alla2b[k][i] == -1) alla2b[k][i] = last_idx + 1;
                else
                    alla2b[k][i] += 1;
            }
//            last_idx += 1;
        }
        last_idx += 1;
        cout << last_idx << endl;
    }
    for(int m=0;m<alla2b.size();m++){
        for(int n=0;n<alla2b[m].size();n++){
            cout << alla2b[m][n] << " ";
        }
        cout << endl;
    }
    cout << "END" << endl;
    for(int t=0;t<alla2b.size();t++){
        last_idx = -1;
        vec_string seqa2m(seqLen+maxLen,"-");

        for(int i=0;i<alla2b[t].size();i++){
//            cout << seqs[t][i] << " ";
            string tem(1,seqs[t][i]);
            if(alla2b0b[t][i] == -1){
                transform(tem.begin(),tem.end(),tem.begin(),::tolower);
                seqa2m[alla2b[t][i]] = tem;
//                last_idx += 1;
            }
            else{
                transform(tem.begin(),tem.end(),tem.begin(),::toupper);
                seqa2m[alla2b[t][i]] = tem;
//                last_idx = alla2b[t][i];
            }
//            cout << tem;
        }
        cout << endl;
        for(int k=0;k<seqa2m.size();k++)
            cout << seqa2m[k];
        cout << endl;
    }
    return alignedSeqs;
}

vec_string updatealnnew1(vector<vec_int> alla2b,vector<string> seqs,string refseq,vec_string allIDs){

    cout << "MSA START:" << endl;
    vector<vec_string> alignedSeqs;
    vec_string alignedSeqs1;
    vector<vec_int> alla2b0b = alla2b;
    int maxLen = 0;
    int seqLen = 0;
    int last_idx = -1;
    for(int m=0;m<alla2b.size();m++){
        if(maxLen < alla2b[m].size()) maxLen = alla2b[m].size();
    }
    int countNum=0;
    for(int i=0;i<maxLen;i++){
        for(int j=0;j<alla2b.size();j++){
            if(i>=alla2b[j].size()) continue;
            if(alla2b[j][i] == -1) {countNum++;break;}
        }
    }
//    for(int m=0;m<alla2b.size();m++){
//        for(int n=0;n<alla2b[m].size();n++){
//            cout << alla2b[m][n] << " ";
//        }
//        cout << endl;
//    }
    for(int i=0;i<maxLen;i++){
        int flag = 0;
        for (int t=0;t<alla2b.size();t++){
            if (i>=alla2b[t].size()) continue;
            if (alla2b[t][i] == -1) {flag=1;break;}
        }
        if(flag){
            if(i==0) {last_idx = -1;}
            else {
                int temidx = 9999;
                for(int t=0;t<alla2b.size();t++){
                    if (i>=alla2b[t].size()) continue;
                    if (alla2b[t][i-1] < temidx) temidx = alla2b[t][i-1];
//                    if(i==24 || i == 25 || i==26) cout << alla2b[t][i-1] << endl;
                }
                last_idx = temidx;
            }
//            cout << i << " : " << last_idx << endl;
            for(int t=0;t<alla2b.size();t++){
                for(int k=0;k<alla2b[t].size();k++){
                    if (alla2b[t][k] > last_idx) alla2b[t][k] += 1;
                }
                if(i>=alla2b[t].size()) continue;
                if(alla2b[t][i] == -1){
                    if(i==0) alla2b[t][i] = last_idx + 1;
                    else alla2b[t][i] = alla2b[t][i-1] + 1;
                }
            }
        }
    }
    /*
    for(int m=0;m<alla2b.size();m++){
        if (alla2b[m].size() > seqLen) seqLen = alla2b[m].size();
        int inslen = 0;
        for(int n=0;n<alla2b[m].size();n++){
            if(alla2b[m][n] == -1) inslen++;
        }
        if(inslen > maxLen) maxLen = inslen;
    }

    for(int m=0;m<alla2b.size();m++){
        for(int n=0;n<alla2b[m].size();n++){
            cout << alla2b[m][n] << " ";
        }
        cout << endl;
    }

    cout << maxLen << endl;
    vec_int tmpList;
    int last_idx = -1;
    int flag = 0;
    for(int i=0;i<seqLen;i++){
        flag = 0;
        for(int j=0;j<alla2b.size();j++){
            cout << alla2b[j][i] << " ";
            if(i>=alla2b[j].size()) continue;
            if(alla2b[j][i] == -1) {flag=1;break;}
        }
        cout << " TEST " << flag;
        if(flag) {
            for(int m=0;m<alla2b.size();m++){
                if (i>=alla2b[m].size()) continue;
                for(int n=i+1;n<alla2b[m].size();n++){
                    if(alla2b[m][n] > last_idx) alla2b[m][n] += 1;
                }
            }
            for(int k=0;k<alla2b.size();k++){
                if (alla2b[k][i] == -1) alla2b[k][i] = last_idx + 1;
                else
                    alla2b[k][i] += 1;
            }
//            last_idx += 1;
        }
        last_idx += 1;
        cout << last_idx << endl;
    }*/
//    for(int m=0;m<alla2b.size();m++){
//        for(int n=0;n<alla2b[m].size();n++){
//            cout << alla2b[m][n] << " ";
//        }
//        cout << endl;
//    }
//    cout << "END" << endl;
//    cout << countNum << ":" << maxLen << endl;
    for(int t=0;t<alla2b.size();t++){
        last_idx = -1;
        vec_string seqa2m(countNum+maxLen,"-");

        for(int i=0;i<alla2b[t].size();i++){
//            cout << seqs[t][i] << " ";
            string tem(1,seqs[t][i]);
            if(alla2b0b[t][i] == -1){
                transform(tem.begin(),tem.end(),tem.begin(),::tolower);
                seqa2m[alla2b[t][i]] = tem;
//                last_idx += 1;
            }
            else{
                transform(tem.begin(),tem.end(),tem.begin(),::toupper);
                seqa2m[alla2b[t][i]] = tem;
//                last_idx = alla2b[t][i];
            }
//            cout << tem;
        }
        alignedSeqs.push_back(seqa2m);
//        cout << endl;
//        for(int k=0;k<seqa2m.size();k++)
//            cout << seqa2m[k];
//        cout << endl;
    }
    vec_int tmpout;
    for(int i=0;i<countNum+maxLen;i++){
        int flag = 0;
        for(int j=0;j<alignedSeqs.size();j++){
            if(alignedSeqs[j][i] != "-") {flag =1;break;}
        }
        if (flag==0) tmpout.push_back(i);
       }
    for(int i=0;i<alignedSeqs.size();i++){
        cout << ">" << allIDs[i] << endl;
        for(int j=0;j<countNum+maxLen;j++){
            if (count(tmpout.begin(),tmpout.end(),j)) continue;
            cout << alignedSeqs[i][j];
        }
        cout << endl;
    }
    return alignedSeqs1;
}

vec_string updatealn(vector<vec_int> alla2b,vector<string> seqs,string refseq){
    int last_idx = -1;
    vector<vec_int> allb2a;
//    cout << "start" << endl;
//    cout << refseq << endl;
    for(int i=0;i<alla2b.size();i++){
        for(int j=0;j<alla2b[i].size();j++){
            cout << alla2b[i][j] << " ";
        }
        cout << endl;
    }
    for(int i=0;i<alla2b.size();i++){
        vec_int tmpb2a(refseq.length(),-1);
        for(int j=0;j<alla2b[i].size();j++){
//            cout << j << ":" << alla2b[i][j] << " ";
            if(alla2b[i][j] == -1) continue;
//            cout << j << ":" << alla2b[i][j] << " ";
            tmpb2a[alla2b[i][j]] = j;
        }
        for(int j=0;j<tmpb2a.size();j++)
            cout << tmpb2a[j] << " ";
        cout << endl;
//            cout << alla2b[i][j] << " ";
//        cout << endl;
        allb2a.push_back(tmpb2a);
    }
    cout << allb2a.size() << endl;
    int inslen = 0;
//    for(int i=0;i<refseq.length();i++){
//        for(int j=0;j<allb2a.size();j++){
////            cout << allb2a[j][i] << " ";
//            if(allb2a[j][i] == -1) {inslen++;break;}
//        }
//    }
//    last_idx = -1;
//    for(int t=0;t<allb2a.size();t++){
//        last_idx = -1;
//        for(int k=0;k<allb2a[t].size();k++){
//            if(allb2a[t][k]==-1) {
//                for(int m=0;m<allb2a.size();m++){
//                    for(int n=0;n<allb2a[m].size();n++){
//                        if(allb2a[m][n] > last_idx) {
//                            allb2a[m][n] += 1;
//                        }
//                    }
//                }
//                allb2a[t][k] = last_idx +1;
//                break;
//            }
//            last_idx = allb2a[t][k];
//        }
//    }
//    for(int t=0;t<allb2a.size();t++){
//        for(int m=0;m<allb2a[t].size();m++)
//            cout << allb2a[t][m] << " ";
//        cout << endl;
//    }
//    for(int t=0;t<alla2b.size();t++){
//        for(int m=0;m<alla2b[t].size();m++)
//            cout << alla2b[t][m] << " ";
//        cout << endl;
//    }
//    for(int i=0;i<refseq.length();i++)
//    for(int i=0;i<refseq.length();i++){
//        for(int j=0;j<allb2a.size();j++){
////            cout << allb2a[j][i] << " ";
//            if(allb2a[j][i] == -1) {inslen++;break;}
//        }
////        cout << endl;
//        for(int j=0;j<allb2a.size();j++){
//            if(allb2a[j][i] == -1) {
//                for(int k1=0;k1<allb2a.size();k1++){
//                    for(int k=i;k<refseq.length();k++) {
//                        if(allb2a[k1][k] > i) allb2a[k1][k] += 1;
//                    }
//                }
//                allb2a[j][i] = last_idx + 1;
//                break;
//            }
////            allb2a[j][i] = last_idx + 1;
//            last_idx = allb2a[j][i];
//        }
////        last_idx =
//    }
//    for(int i=0;i<refseq.length();i++)
    vector<vec_int> alla2b0b = alla2b;
//    cout << "inslen:" << inslen << endl;
    vec_string alignedSeqs;
    int maxLen=0;
//    for(int m=0;m<alla2b.size();m++){
//        last_idx = -1;
//        inslen = 0;
//        for(int i=0;i<alla2b[m].size();i++){
//            if (alla2b[m][i] == -1){
//                for(int j=i+1;j<alla2b[m].size();j++){
//                    if(alla2b[m][j]>last_idx)
//                        alla2b[m][j] += 1;
//                }
//                alla2b[m][i] = last_idx + 1;
//                inslen++;
//            }
//            last_idx = alla2b[m][i];
//        }
//        if (inslen > maxLen) maxLen = inslen;
//    }


     for(int m=0;m<alla2b.size();m++){
        last_idx = -1;
        inslen = 0;
        for(int i=0;i<alla2b[m].size();i++){
            if (alla2b[m][i] == -1){
                for(int k=0;k<alla2b.size();k++){
                    for(int j=i+1;j<alla2b[m].size();j++){
                        if(alla2b[m][j]>last_idx)
                            alla2b[m][j] += 1;
                    }
                    alla2b[m][i] = last_idx + 1;
                    inslen++;
                }
            }
            last_idx = alla2b[m][i];
        }
        if (inslen > maxLen) maxLen = inslen;
    }



    for(int m=0;m<alla2b.size();m++){
        for(int n=0;n<alla2b[m].size();n++){
            cout << alla2b[m][n] << " ";
        }
        cout << endl;
    }

//    int last_idx = -1;
    for(int t=0;t<alla2b.size();t++){
        last_idx = -1;
        vec_string seqa2m(refseq.length()+maxLen,"-");

        for(int i=0;i<alla2b[t].size();i++){
//            cout << seqs[t][i] << " ";
            string tem(1,seqs[t][i]);
            if(alla2b0b[t][i] == -1){
                transform(tem.begin(),tem.end(),tem.begin(),::tolower);
                seqa2m[alla2b[t][i]] = tem;
//                last_idx += 1;
            }
            else{
                transform(tem.begin(),tem.end(),tem.begin(),::toupper);
                seqa2m[alla2b[t][i]] = tem;
//                last_idx = alla2b[t][i];
            }
//            cout << tem;
        }
        cout << endl;
        for(int k=0;k<seqa2m.size();k++)
            cout << seqa2m[k];
        cout << endl;
    }

//    for(int t=0;t<alla2b.size();t++){
//        last_idx = -1;
//        for(int i=0;i<alla2b[t].size();i++){
//            if(alla2b[t][i] == -1)
//                for(int j=0;j<alla2b.size();j++)
//                    for(int l=0;l<alla2b[j].size();l++){
//                        if(alla2b[j][l] > last_idx)
//                            alla2b[j][l] += 1;
//                }
//            last_idx = alla2b[t][i];
//        }
//    }
//    int maxLen = 0;
//    cout << "maxLen" << maxLen << endl;
//    for(int t=0;t<alla2b.size();t++)
//        for(int i=0;i<alla2b[t].size();i++){
//            alla2b[t][i] += 1;
//            if(alla2b[t][i] > maxLen) maxLen = alla2b[t][i];
//        }
//    vec_string alignedSeqs;
//    for(int i=0;i<alla2b.size();i++){
//        string temseq(maxLen,'-');
//        for(int j=0;j<alla2b[i].size();j++){
//            cout << alla2b[i][j] << " ";
//            if(alla2b[i][j] > 0)
////                temseq += seqs[i][alla2b[i][j]];
//                temseq += seqs[i][j];
//            else
//                temseq += "-";
//        }
//        cout << endl;
//        alignedSeqs.push_back(temseq);
//    }
    return alignedSeqs;
}

void pair_a2b2aln(vec_int a2b1,vec_int a2b2,string seqs){
    int last_idx = -1;
    for(int i=0;i<a2b1.size();i++){
        if(a2b1[i] == -1){
            for(int j=0;j<a2b1.size();j++){
                if (a2b1[i] > last_idx){
                    a2b1[j] += 1;
                    a2b2[j] += 1;
                }
            }
        }
        last_idx = a2b1[i];
    }
    last_idx = -1;
    for(int i=0;i<a2b2.size();i++){
        if(a2b2[i] == -1){
            for(int j=0;j<a2b2.size();j++){
                if (a2b2[i] > last_idx){
                    a2b1[j] += 1;
                    a2b2[j] += 1;
                }
            }
        }
        last_idx = a2b2[i];
    }

}
