/*************************************************************************
    > File Name: msa.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Wed May 25 09:22:16 2022
 ************************************************************************/

#include<iostream>
#include "msa.h"
using namespace std;

int readMSAsto(string msaFile,vector<string> &Names, vector<string> &seqAlns, string &RNAtype, string &RNAss){
    // read stockholm file
    const char *aliFile;
    aliFile = msaFile.c_str();
    char ** names;
    char ** seqalns;
    char * id;
    char * structures;
    int alnnum;
    alnnum = vrna_file_msa_read(aliFile,&names,&seqalns,&id,&structures,VRNA_FILE_FORMAT_MSA_STOCKHOLM);
    if(alnnum == 0) {
        cout << "stockholm msa file is empty!" << endl;
        return 0;
    }
    int i = 0;
    while(1){
        if(seqalns[i]){
//            cout << names[i] << endl;
//            cout << seqalns[i] << endl;
          string nametmp = names[i];
          Names.push_back(nametmp);
          string seqtmp = seqalns[i];
          seqAlns.push_back(seqtmp);
        }else break;
        i++;
    }
//	cout << structures << endl;
	if (structures == NULL) RNAss = "UnknownSS";
	else RNAss = structures;
	if (id == NULL) RNAtype = "Unknown";
	else RNAtype = id;
	cout << RNAss << " " << RNAtype << endl;
//    RNAtype = id;
    return alnnum;
}

vec_string msafam2prf(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns){
    int seqNum=0;
    string seqtmp;
    cout << "start with function" <<endl;
    for(int i=0;i<msanames.size();i++){
        if(!queryName.compare(msanames[i])) {
            //string temseq = replace_all(msaalns[i],'-','');
            string temseq = msaalns[i];
            seqtmp = msaalns[i];
			cout << temseq << endl;
            temseq.erase(remove(temseq.begin(),temseq.end(),'-'),temseq.end());
            if(temseq.length() != querySeq.length()) {cout << "sequence is not equal!!";exit(1);}
            int k =0;
            for(int j=0;j<msaalns[i].size();j++){
                if(msaalns[i][j] != '-' && msaalns[i][j] == querySeq[k] ) k++;
            }
            if(querySeq.length() != k) {cout << "sequence is not equal!!!";exit(1);}
            seqNum = i;
            break;
        }
    }
    cout << "start with prfM" << endl;
    vector<vector<float>> prfM(querySeq.length(),vector<float>(4,0));
    cout << "flag" << endl;
    cout << "seqNum " << seqNum << " querySeqlen" << querySeq.length()<< endl;
    cout << prfM.size() << " " << prfM[0].size() << endl;
    int t = 0;
//  string temseq = msaalns[seqNum];
    cout << seqtmp << endl;
    cout << "TEST" << endl;
    cout << seqtmp.length() << endl;
    for(int i=0;i<seqtmp.length();i++){
        int cA=0,cU=0,cC=0,cG=0;
        cout << msaalns[seqNum][i] << " " ;
        if (msaalns[seqNum][i] == '-') continue;
        for(int k =0;k<msaalns.size();k++){ 
            if(msaalns[k][i] == 'A') cA++;
            else if(msaalns[k][i] == 'U') cU++;
            else if(msaalns[k][i] == 'C') cC++;
            else if(msaalns[k][i] == 'G') cG++;
            else continue;
        }
//        cout << "t " << t << " seqlen " << querySeq.length() << endl;
        prfM[t][0] = 1.0 * cA / (cA + cU + cC + cG);
        prfM[t][1] = 1.0 * cU / (cA + cU + cC + cG);
        prfM[t][2] = 1.0 * cC / (cA + cU + cC + cG);
        prfM[t][3] = 1.0 * cG / (cA + cU + cC + cG);
        t++;
    }
    cout << "start with prflists" << endl;
    vec_string prflists;
    int flag=0;
    for(int i=0;i<seqtmp.length();i++){
        if(msaalns[seqNum][i] == '-') continue;
        if(msaalns[seqNum][i] == 'N') continue;
        ostringstream os1;
        if(RNAss[i] == '.')
            os1 << "PRF " << flag << " " << querySeq[flag] << " N\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        else
            os1 << "PRF " << flag << " " << querySeq[flag] << " P\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        prflists.push_back(os1.str());
        os1.clear();
        flag ++;
    }
    return prflists;

}

vec_string msafam2prf(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns,string type){
    int seqNum=0;
    string seqtmp = querySeq;
    cout << "start with function" <<endl;
//    for(int i=0;i<msanames.size();i++){
//        if(!queryName.compare(msanames[i])) {
//            //string temseq = replace_all(msaalns[i],'-','');
//            string temseq = msaalns[i];
//            seqtmp = msaalns[i];
//			cout << temseq << endl;
//            temseq.erase(remove(temseq.begin(),temseq.end(),'-'),temseq.end());
//            if(temseq.length() != querySeq.length()) {cout << "sequence is not equal!!";exit(1);}
//            int k =0;
//            for(int j=0;j<msaalns[i].size();j++){
//                if(msaalns[i][j] != '-' && msaalns[i][j] == querySeq[k] ) k++;
//            }
//            if(querySeq.length() != k) {cout << "sequence is not equal!!!";exit(1);}
//            seqNum = i;
//            break;
//        }
//    }
//    cout << "start with prfM" << endl;
    vector<vector<float>> prfM(querySeq.length(),vector<float>(4,0));
    cout << "flag" << endl;
    cout << "seqNum " << seqNum << " querySeqlen" << querySeq.length()<< endl;
    cout << prfM.size() << " " << prfM[0].size() << endl;
    int t = 0;
//  string temseq = msaalns[seqNum];
//    cout << seqtmp << endl;
//    cout << "TEST" << endl;
//    cout << seqtmp.length() << endl;
    for(int i=0;i<seqtmp.length();i++){
        int cA=0,cU=0,cC=0,cG=0;
//        cout << msaalns[seqNum][i] << " " ;
//        if (msaalns[seqNum][i] == '-') continue;
        for(int k =0;k<msaalns.size();k++){
            if(msaalns[k][i] == 'A') cA++;
            else if(msaalns[k][i] == 'U') cU++;
            else if(msaalns[k][i] == 'C') cC++;
            else if(msaalns[k][i] == 'G') cG++;
            else continue;
        }
        if(cA+cU+cC+cG == 0) {
            prfM[t][0] = 0.25;
            prfM[t][1] = 0.25;
            prfM[t][2] = 0.25;
            prfM[t][3] = 0.25;
            t++;
            continue;
        }
//        cout << "t " << t << " seqlen " << querySeq.length() << endl;
        prfM[t][0] = 1.0 * cA / (cA + cU + cC + cG);
        prfM[t][1] = 1.0 * cU / (cA + cU + cC + cG);
        prfM[t][2] = 1.0 * cC / (cA + cU + cC + cG);
        prfM[t][3] = 1.0 * cG / (cA + cU + cC + cG);
        t++;
    }
    cout << "start with prflists" << endl;
    vec_string prflists;
    int flag=0;
    for(int i=0;i<seqtmp.length();i++){
//        if(msaalns[seqNum][i] == '-') continue;
//        if(msaalns[seqNum][i] == 'N') continue;
//        if(seqtmp[i] == '-') continue;
//        if(seqtmp[i] == 'N') continue;
        ostringstream os1;
        if(seqtmp[i] == '.' || seqtmp[i] == 'N')
            os1 << "PRF " << flag << " " << querySeq[flag] << " N\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        else {
            if(RNAss[i] == '.')

                os1 << "PRF " << flag << " " << querySeq[flag] << " N\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
            else
                os1 << "PRF " << flag << " " << querySeq[flag] << " P\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
             }
        prflists.push_back(os1.str());
        os1.clear();
        flag ++;
    }
    return prflists;

}

vec_string msa2prf(string querySeq,string msaFile,string ss){
    vector<string> alnNames;
    vector<string> alnSeqs;
    string RNAtype;
    string RNAss;
    int alnNum = readMSAsto(msaFile,alnNames,alnSeqs,RNAtype,RNAss);
    if(alnNum==0){
        cout << "MSA file is empty!" << endl;
    }
    int j = 0;
    for(int i=0;i<alnSeqs[0].size();i++){
        if(alnSeqs[0][i] != '-' && alnSeqs[0][i] == querySeq[j]) {
            j++;
        }
        else{
            cout << "QuerySeq is not same as the sequence in MSA results!" << endl;
        }
    }
    vector<vector<float>> prfM(querySeq.length(),vector<float>(4,0));
    int t = 0;
    for(int i=0;i<alnSeqs[0].size();i++){
        int cA=0,cU=0,cC=0,cG=0;
        if (alnSeqs[0][t] == '-') continue;
        for(int k =0;k<alnNum;k++){
            if(alnSeqs[k][i] == 'A') cA += 1;
            else if(alnSeqs[k][i] == 'U') cU += 1;
            else if(alnSeqs[k][i] == 'C') cC += 1;
            else if(alnSeqs[k][i] == 'G') cG += 1;
            else continue;
        }
        prfM[t][0] = 1.0 * cA / (cA + cU + cC + cG);
        prfM[t][1] = 1.0 * cU / (cA + cU + cC + cG);
        prfM[t][2] = 1.0 * cC / (cA + cU + cC + cG);
        prfM[t][3] = 1.0 * cG / (cA + cU + cC + cG);
        t++;
    }
    vec_string prflists;
    int flag=0;
    for(int i=0;i<alnSeqs[0].size();i++){
        if(alnSeqs[0][i] == '-') continue;
        if(alnSeqs[0][i] == 'N') continue;
        ostringstream os1;
        if(RNAss[i] == '.')
            os1 << "PRF " << flag << " " << querySeq[flag] << " N\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        else
            os1 << "PRF " << flag << " " << querySeq[flag] << " P\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        prflists.push_back(os1.str());
        os1.clear();
        flag ++;
    }
    return prflists;
}

vec_string msafam2prfWithWeight(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns,float weight){
    int seqNum=0;
    int alnNum = msaalns.size();
    string seqtmp;
    cout << "start with function" <<endl;
    for(int i=0;i<msanames.size();i++){
        if(!queryName.compare(msanames[i])) {
            //string temseq = replace_all(msaalns[i],'-','');
            string temseq = msaalns[i];
            seqtmp = msaalns[i];
            temseq.erase(remove(temseq.begin(),temseq.end(),'-'),temseq.end());
            if(temseq.length() != querySeq.length()) {cout << "sequence is not equal!!";exit(1);}
            int k =0;
            for(int j=0;j<msaalns[i].size();j++){
                if(msaalns[i][j] != '-' && msaalns[i][j] == querySeq[k] ) k++;
            }
            if(querySeq.length() != k) {cout << "sequence is not equal!!!";exit(1);}
            seqNum = i;
            break;
        }
    }

    Eigen::MatrixXf Neff = Eigen::MatrixXf::Zero(alnNum,1);
    Eigen::MatrixXf Init_prf = Eigen::MatrixXf::Zero(querySeq.length(),4);
    for(int i=0;i<alnNum;i++)
        for(int j=0;j<alnNum;j++){
            if( HammingDistance(msaalns[i],msaalns[j]) <= 0.2 ) Neff(i,0) = Neff(i,0) + 1;
        }
    for(int i=0;i<alnNum;i++)
        Neff(i,0) = 1.0 / Neff(i,0);

    float prf_tmp = 0.6;
    float prf_tmp_o = (1-prf_tmp) / 3;
    for(int i=0;i<querySeq.length();i++){
        for(int j=0;j<4;j++){
            Init_prf(i,j) = prf_tmp_o;
        }
        if(querySeq[i] == 'A') Init_prf(i,0) = prf_tmp;
        else if(querySeq[i] == 'U') Init_prf(i,1) = prf_tmp;
        else if(querySeq[i] == 'C') Init_prf(i,2) = prf_tmp;
        else if(querySeq[i] == 'G') Init_prf(i,3) = prf_tmp;
        else continue;
    }

    cout << "start with prfM" << endl;
    vector<vector<float>> prfM(querySeq.length(),vector<float>(4,0));
    cout << "flag" << endl;
    cout << "seqNum " << seqNum << " querySeqlen" << querySeq.length()<< endl;
    cout << prfM.size() << " " << prfM[0].size() << endl;
    int t = 0;
//  string temseq = msaalns[seqNum];
    cout << seqtmp << endl;
    cout << "TEST" << endl;
    cout << seqtmp.length() << endl;
    for(int i=0;i<seqtmp.length();i++){
        int cA=0,cU=0,cC=0,cG=0;
        cout << msaalns[seqNum][i] << " " ;
        if (msaalns[seqNum][i] == '-') continue;
        for(int k =0;k<msaalns.size();k++){ 
            if(msaalns[k][i] == 'A') cA += Neff(k,0);
            else if(msaalns[k][i] == 'U') cU += Neff(k,0);
            else if(msaalns[k][i] == 'C') cC += Neff(k,0);
            else if(msaalns[k][i] == 'G') cG += Neff(k,0);
            else continue;
        }
//        cout << "t " << t << " seqlen " << querySeq.length() << endl;
        prfM[t][0] = 1.0 * cA / (cA + cU + cC + cG) * weight + (1-weight) * Init_prf(t,0);
        prfM[t][1] = 1.0 * cU / (cA + cU + cC + cG) * weight + (1-weight) * Init_prf(t,1);
        prfM[t][2] = 1.0 * cC / (cA + cU + cC + cG) * weight + (1-weight) * Init_prf(t,2);
        prfM[t][3] = 1.0 * cG / (cA + cU + cC + cG) * weight + (1-weight) * Init_prf(t,3);
        t++;
    }
    cout << "start with prflists" << endl;
    vec_string prflists;
    int flag=0;
    for(int i=0;i<seqtmp.length();i++){
        if(msaalns[seqNum][i] == '-') continue;
        if(msaalns[seqNum][i] == 'N') continue;
        ostringstream os1;
        if(RNAss[i] == '.')
            os1 << "PRF " << flag << " " << querySeq[flag] << " N\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        else
            os1 << "PRF " << flag << " " << querySeq[flag] << " P\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        prflists.push_back(os1.str());
        os1.clear();
        flag ++;
    }
    return prflists;

}

vec_string msafam2prfWithWeight(string querySeq,string queryName,string RNAss,vec_string msanames,vec_string msaalns){
    int seqNum=0;
    int alnNum = msaalns.size();
    string seqtmp;
    cout << "start with function" <<endl;
    for(int i=0;i<msanames.size();i++){
        if(!queryName.compare(msanames[i])) {
            //string temseq = replace_all(msaalns[i],'-','');
            string temseq = msaalns[i];
            seqtmp = msaalns[i];
            temseq.erase(remove(temseq.begin(),temseq.end(),'-'),temseq.end());
            if(temseq.length() != querySeq.length()) {cout << "sequence is not equal!!";exit(1);}
            int k =0;
            for(int j=0;j<msaalns[i].size();j++){
                if(msaalns[i][j] != '-' && msaalns[i][j] == querySeq[k] ) k++;
            }
            if(querySeq.length() != k) {cout << "sequence is not equal!!!";exit(1);}
            seqNum = i;
            break;
        }
    }

    Eigen::MatrixXf Neff = Eigen::MatrixXf::Zero(alnNum,1);
    for(int i=0;i<alnNum;i++)
        for(int j=0;j<alnNum;j++){
            if( HammingDistance(msaalns[i],msaalns[j]) <= 0.2 ) Neff(i,0) = Neff(i,0) + 1;
        }
    for(int i=0;i<alnNum;i++)
        Neff(i,0) = 1.0 / Neff(i,0);

    cout << "start with prfM" << endl;
    vector<vector<float>> prfM(querySeq.length(),vector<float>(4,0));
    cout << "flag" << endl;
    cout << "seqNum " << seqNum << " querySeqlen" << querySeq.length()<< endl;
    cout << prfM.size() << " " << prfM[0].size() << endl;
    int t = 0;
//  string temseq = msaalns[seqNum];
    cout << seqtmp << endl;
    cout << "TEST" << endl;
    cout << seqtmp.length() << endl;
    for(int i=0;i<seqtmp.length();i++){
        int cA=0,cU=0,cC=0,cG=0;
        cout << msaalns[seqNum][i] << " " ;
        if (msaalns[seqNum][i] == '-') continue;
        for(int k =0;k<msaalns.size();k++){ 
            if(msaalns[k][i] == 'A') cA += Neff(k,0);
            else if(msaalns[k][i] == 'U') cU += Neff(k,0);
            else if(msaalns[k][i] == 'C') cC += Neff(k,0);
            else if(msaalns[k][i] == 'G') cG += Neff(k,0);
            else continue;
        }
//        cout << "t " << t << " seqlen " << querySeq.length() << endl;
        prfM[t][0] = 1.0 * cA / (cA + cU + cC + cG);
        prfM[t][1] = 1.0 * cU / (cA + cU + cC + cG);
        prfM[t][2] = 1.0 * cC / (cA + cU + cC + cG);
        prfM[t][3] = 1.0 * cG / (cA + cU + cC + cG);
        t++;
    }
    cout << "start with prflists" << endl;
    vec_string prflists;
    int flag=0;
    for(int i=0;i<seqtmp.length();i++){
        if(msaalns[seqNum][i] == '-') continue;
        if(msaalns[seqNum][i] == 'N') continue;
        ostringstream os1;
        if(RNAss[i] == '.')
            os1 << "PRF " << flag << " " << querySeq[flag] << " N\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        else
            os1 << "PRF " << flag << " " << querySeq[flag] << " P\t" << prfM[flag][0] << "\t" <<prfM[flag][1] << "\t" << prfM[flag][2]  << "\t" << prfM[flag][3] << endl;
        prflists.push_back(os1.str());
        os1.clear();
        flag ++;
    }
    return prflists;

}
