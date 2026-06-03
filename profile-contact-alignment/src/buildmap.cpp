/*************************************************************************
    > File Name: buildmap.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun May 15 16:20:15 2022
 ************************************************************************/

#include "buildmap.h"
//#include "msa.h"
#include<iostream>
using namespace std;


Matrix4d PAM(int num){
    Matrix4d results;
    Matrix4d initMat;
    results << 0.99,0.0033,0.0033,0.0033,
               0.0033,0.99,0.0033,0.0033,
               0.0033,0.0033,0.99,0.0033,
               0.0033,0.0033,0.0033,0.99;
    initMat << 0.99,0.0033,0.0033,0.0033,
               0.0033,0.99,0.0033,0.0033,
               0.0033,0.0033,0.99,0.0033,
               0.0033,0.0033,0.0033,0.99;
    for(int i=1;i<num;i++)
        results *= initMat;
    return results;
}

vec_string seq2prf(string seq,float prob,string ss){
    // sequence to prf
    vec_string prflists;
    float basePro = (1.0 - prob) / 3;
    
    for(int i=0;i<seq.length();i++)
    {
        float pA = basePro;
        float pU = basePro;
        float pC = basePro;
        float pG = basePro;
        if(seq[i]=='A')
            pA = prob;
        else if (seq[i] == 'U')
            pU = prob;
        else if (seq[i] == 'C')
            pC = prob;
        else if (seq[i] == 'G')
            pG = prob;
        else
            continue;
        ostringstream os1;
        if (ss[i] == '.')
            os1 <<"PRF "<<i << " " << seq[i] << " N\t" << pA << "\t" << pU << "\t" << pC << "\t" << pG << endl;
        else
            os1 <<"PRF "<<i << " " << seq[i] << " P\t" << pA << "\t" << pU << "\t" << pC << "\t" << pG << endl;
        prflists.push_back(os1.str());
        os1.clear();
    }
    return prflists;
}

vec_string seq2prf(string seq, float prob){
    // sequence to prf
    vec_string prflists;
    float basePro = (1.0 - prob) / 3;
    
    for(int i=0;i<seq.length();i++)
    {
        float pA = basePro;
        float pU = basePro;
        float pC = basePro;
        float pG = basePro;
        if(seq[i]=='A')
            pA = prob;
        else if (seq[i] == 'U')
            pU = prob;
        else if (seq[i] == 'C')
            pC = prob;
        else if (seq[i] == 'G')
            pG = prob;
        else
            continue;
        ostringstream os1;
        os1 <<"PRF "<<i << " " << seq[i] << " N\t" << pA << "\t" << pU << "\t" << pC << "\t" << pG << endl;
        prflists.push_back(os1.str());
        os1.clear();
    }
    return prflists;
}

vec_string seq2map(string seq){
    // sequence to contact map
    vec_string maps;
    ostringstream os1;
    os1 << "LEN " << seq.length() << endl;
//    cout << "seq.length()" << seq.length() << endl;
    maps.push_back(os1.str());
    os1.clear();
    for(int i=0;i<seq.length();i++)
        for (int j=i+1;j<seq.length();j++)
        {
            ostringstream os1;
            if((seq[i] == 'A' && seq[j] == 'U') || ( seq[i] == 'U' && seq[j] == 'A') ||  (seq[i] == 'C' && seq[j] == 'G') ||  (seq[i] == 'G' && seq[j] == 'C'))
                os1 << "CON " << i << " " << j << " 1.0" << endl;
            else if((seq[i] == 'G' && seq[j] == 'U') || ( seq[i] == 'U' && seq[j] == 'G'))
                os1 << "CON " << i << " " << j << " 0.1" << endl;
            if (os1.str().empty()) continue;
            maps.push_back(os1.str());
            os1.clear();
        }
    return maps;
}

vec_string ss2map(string ss){
    // secondary structure to contact map 
    vec_string ssmap;
    ostringstream os1;
    os1 <<"LEN " << ss.length() << endl;
    ssmap.push_back(os1.str());
    os1.clear();
    stack<int> ss1_l,ss2_l,ss3_l,ss4_l;
    stack<int> ss5_l,ss6_l,ss7_l,ss8_l;

    cout <<"ss length is " << ss.length() << endl;
    for(int i=0;i<ss.length();i++){
        if(ss[i] == '(')
            ss1_l.push(i);
        else if(ss[i] == '<')
            ss2_l.push(i);
        else if(ss[i] == '[')
            ss3_l.push(i);
        else if(ss[i] == '{')
            ss4_l.push(i);
        else if(ss[i] == 'A')
            ss5_l.push(i);
        else if(ss[i] == 'B')
            ss6_l.push(i);
        else if(ss[i] == 'C')
            ss7_l.push(i);
        else if(ss[i] == 'D')
            ss8_l.push(i);
        else if(ss[i] == 'a'){
            ostringstream os1;
            os1 << "CON " << ss5_l.top() << " " << i << " 1.0"<< endl;
            ss5_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else if(ss[i] == 'b'){
            ostringstream os1;
            os1 << "CON " << ss6_l.top() << " " << i << " 1.0"<< endl;
            ss6_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else if(ss[i] == 'c'){
            ostringstream os1;
            os1 << "CON " << ss7_l.top() << " " << i << " 1.0"<< endl;
            ss7_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else if(ss[i] == 'd'){
            ostringstream os1;
            os1 << "CON " << ss8_l.top() << " " << i << " 1.0"<< endl;
            ss8_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else if(ss[i] == ')'){
            ostringstream os1;
            os1 << "CON " << ss1_l.top() << " " << i << " 1.0"<< endl;
            ss1_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else if(ss[i] == '>'){
            ostringstream os1;
            os1 << "CON " << ss2_l.top() << " " << i << " 1.0"<< endl;
            ss2_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else if(ss[i] == ']'){
            ostringstream os1;
            os1 << "CON " << ss3_l.top() << " " << i << " 1.0"<< endl;
            ss3_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();

        }
        else if(ss[i] == '}'){
            ostringstream os1;
            os1 << "CON " << ss4_l.top() << " " << i << " 1.0"<< endl;
            ss4_l.pop();
            ssmap.push_back(os1.str());
            os1.clear();
        }
        else
            continue;
    }
    if (ss1_l.size() != 0 or ss2_l.size() != 0 or ss3_l.size() != 0 or ss4_l.size() != 0)
    {
        cout << "SS have an error!!!" << endl;
        exit(1);
    }
    std::sort(ssmap.begin()+1,ssmap.end());
    return ssmap;
}

void alignseqs(string seq1,string seq2, vector<vector<int>> maps){

    string news1,news2;
    string s1start,s2start;
    string s1end,s2end;
    int s1startIndex = maps[0][0];
    int s2startIndex = maps[1][0];
    int s1endIndex,s2endIndex;

    int old_i1,old_i2;

    for (int i=0 ; i<maps[0].size();i++)
    {
        if (i==maps[0].size()-1){
            s1endIndex = maps[0][i];
            s2endIndex = maps[1][i];
        }
        if (i){
            for(int j = old_i2+1; j < maps[1][i];j++)
            {
                news1 += "-";
                news2 += seq2[j];
            }
            for(int j = old_i1+1; j < maps[0][i];j++)
            {
                news1 += seq1[j];
                news2 += "-";
            }
        }
        news1 += seq1[maps[0][i]];
        news2 += seq2[maps[1][i]];
        old_i1 = maps[0][i];
        old_i2 = maps[1][i];
    }
    for (int i=0;i<s2startIndex;i++){
        if (i == s2startIndex){
            continue;
        }
        s1start += "-";
        s2start += seq2[i];
    }
    for (int i=0;i<s1startIndex;i++){
        if (i == s1startIndex) continue;
        s1start += seq1[i];
        s2start += "-";
    }
    
    for (int i=s1endIndex+1;i<seq1.length();i++){
        if(i > seq1.length()-1) break;
        s1end += seq1[i];
        s2end += "-";
        if(i == seq1.length()-1) break;
    }
    for (int i=s2endIndex+1;i<seq2.length();i++){
        if (i > seq2.length()-1) break;
        s1end += "-";
        s2end += seq2[i];
        if (i == seq2.length()-1) break;
    }
    news1 = s1start + news1 + s1end;
    news2 = s2start + news2 + s2end;
    cout << "seqA:" << news1 << endl;
    cout << "seqB:" << news2 << endl;
//    cout << ">Identity=" << 1.0 * maps[0].size() / news1.length() << " Length=" << news1.length() << endl;
}

vec_int Map::buildmaps(string &seq,string &seqss, int sep_cutoff,float prob){
    vec_int n2m;
    vec_int m2n;
//    vec_string os_map;
//    vec_string os_prf;
    os_map = ss2map(seqss);
    os_prf = seq2prf(seq,prob,seqss);
    vec_bool range;
    for (int k=0;k<os_map.size();k++)
    {
        stringstream is(os_map[k]);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE")
        {
            int size; is >> size;
            //if(range.size() > 0) {
            //  range.resize(size,0);
            //}
            //else range.resize(size,1);
            range.resize(size,1);

            int m = 0;n2m.resize(size,-1);
            for(int n=0;n<size;n++){
                if(range[n] == 1)
                {
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0));
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else
        {
            cout << os_map[k] << endl;
            cout << label << endl;
            cout << "Contact error!!!" << endl;
            exit(1);
        }

    }
    for (int k =0 ;k <os_prf.size();k++)
    {
        istringstream is(os_prf[k]);
        string label;
        is >> label;
        if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
        else{
            cout << "PRF error!!!" << endl;
            exit(1);
        }

    }
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}

vec_int Map::buildmapswithMSA(string &seq,string &seqName, string &seqss,vec_string msanames,vec_string msaalns, int sep_cutoff){
    vec_int n2m;
    vec_int m2n;
//    vec_string os_map;
//    vec_string os_prf;
    os_map = ss2map(seqss);
    cout << "start with msafam2prf" << endl;
    os_prf = msafam2prf(seq,seqName,seqss,msanames,msaalns,"tmp");
    cout << "end with msafam2prf" << endl;
    vec_bool range;
    for (int k=0;k<os_map.size();k++)
    {
        stringstream is(os_map[k]);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE")
        {
            int size; is >> size;
            //if(range.size() > 0) {
            //  range.resize(size,0);
            //}
            //else range.resize(size,1);
            range.resize(size,1);

            int m = 0;n2m.resize(size,-1);
            for(int n=0;n<size;n++){
                if(range[n] == 1)
                {
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0));
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else
        {
            cout << os_map[k] << endl;
            cout << label << endl;
            cout << "Contact error!!!" << endl;
            exit(1);
        }

    }
    for (int k =0 ;k <os_prf.size();k++)
    {
        istringstream is(os_prf[k]);
        string label;
        is >> label;
        if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
        else{
            cout << "PRF error!!!" << endl;
            exit(1);
        }

    }
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}

vec_int Map::buildmapsMSA(string &seq,string &seqName, string &seqss,vec_string msanames,vec_string msaalns, int sep_cutoff){
    vec_int n2m;
    vec_int m2n;
//    vec_string os_map;
//    vec_string os_prf;
    cout << seq << endl;
    os_map = ss2map(seqss);
    cout << "start with msafam2prf" << endl;
    os_prf = msafam2prf(seq,seqName,seqss,msanames,msaalns);
    cout << "end with msafam2prf" << endl;
    vec_bool range;
    for (int k=0;k<os_map.size();k++)
    {
        stringstream is(os_map[k]);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE")
        {
            int size; is >> size;
            //if(range.size() > 0) {
            //  range.resize(size,0);
            //}
            //else range.resize(size,1);
            range.resize(size,1);

            int m = 0;n2m.resize(size,-1);
            for(int n=0;n<size;n++){
                if(range[n] == 1)
                {
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0));
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else
        {
            cout << os_map[k] << endl;
            cout << label << endl;
            cout << "Contact error!!!" << endl;
            exit(1);
        }

    }
    for (int k =0 ;k <os_prf.size();k++)
    {
        istringstream is(os_prf[k]);
        string label;
        is >> label;
        if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
        else{
            cout << "PRF error!!!" << endl;
            exit(1);
        }

    }
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}

vec_int Map::buildmapsMSA(string &seq,string &seqName, vec_string msanames,vec_string msaalns, int sep_cutoff){
    vec_int n2m;
    vec_int m2n;
//    vec_string os_map;
//    vec_string os_prf;
    string seqss(seq.length(),'-');

    os_map = seq2map(seq);
    cout << "start with msafam2prf" << endl;
    os_prf = msafam2prf(seq,seqName,seqss,msanames,msaalns);
//    os_prf = msafam2prf(seq,seqName,seqss,msanames,msaalns,"tmp");
    cout << "end with msafam2prf" << endl;
    vec_bool range;
    for (int k=0;k<os_map.size();k++)
    {
        stringstream is(os_map[k]);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE")
        {
            int size; is >> size;
            //if(range.size() > 0) {
            //  range.resize(size,0);
            //}
            //else range.resize(size,1);
            range.resize(size,1);

            int m = 0;n2m.resize(size,-1);
            for(int n=0;n<size;n++){
                if(range[n] == 1)
                {
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0));
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else
        {
            cout << os_map[k] << endl;
            cout << label << endl;
            cout << "Contact error!!!" << endl;
            exit(1);
        }

    }
    for (int k =0 ;k <os_prf.size();k++)
    {
        cout << os_prf[k] << endl;
        istringstream is(os_prf[k]);
        string label;
        is >> label;
        if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
        else{
            cout << "PRF error!!!" << endl;
            exit(1);
        }

    }
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}


vec_int Map::load_data_from_file (string file, string &seq, int sep_cutoff,float prob, string type){
    vec_int n2m;
    vec_int m2n;
//    vec_string os_map;
//    vec_string os_prf;
    vec_bool range;
    if(type == "fasta")
    {
        sequence *loadseq = new sequence;
        loadseq = readsinglefasta(file);
        cout << "Database:" << loadseq->id << endl;
        os_map = seq2map(loadseq->seq);
        os_prf = seq2prf(loadseq->seq,prob);
        seq = loadseq->seq;
    }
    else if(type == "dbn")
    {
        seqdbn *queryseqss;
        queryseqss = readdbnfile(file);
        os_map = ss2map(queryseqss->ss);
        //os_prf = seq2prf(queryseqss->seq,prob);
        os_prf = seq2prf(queryseqss->seq,prob,queryseqss->ss);
        //cout << "dbn seq " << queryseqss->seq << endl;
        //cout << "dbn ss " << queryseqss->ss << endl;
        cout << "Query:" << queryseqss->id << endl;
        seq = queryseqss->seq;
    }
    for (int k=0;k<os_map.size();k++)
    {
        stringstream is(os_map[k]);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE")
        {
            int size; is >> size;
            //if(range.size() > 0) {
            //  range.resize(size,0);
            //}
            //else range.resize(size,1);
            range.resize(size,1);
            int m = 0;n2m.resize(size,-1);
            for(int n=0;n<size;n++){
                if(range[n] == 1)
                {
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0));
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else
        {
            cout << os_map[k] << endl;
            cout << label << endl;
            cout << "Contact error!!!" << endl;
            exit(1);
        }

    }
    for (int k =0 ;k <os_prf.size();k++)
    {
        istringstream is(os_prf[k]);
        string label;
        is >> label;
        if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
        else{
            cout << "PRF error!!!" << endl;
            exit(1);
        }
    }
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}

//vec_int Map::seq2data(string seq,int sep_cutoff,vec_bool &range,float prob,string type){
vec_int Map::seq2data(string seq,int sep_cutoff,float prob,string type){
    vec_int n2m;
    vec_int m2n;
//    vec_string os_map;
//    vec_string os_prf;
    
    os_map = seq2map(seq);
    os_prf = seq2prf(seq,prob);
    vec_bool range;

    for (int k=0;k<os_map.size();k++)
    {
        stringstream is(os_map[k]);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE")
        {
            int size; is >> size;
            //cout << size << "size " << range.size() <<endl;
            //if(range.size() > 0) {
            //  range.resize(size,0);
            //}
            //else range.resize(size,1);
            range.resize(size,1);

            int m = 0;n2m.resize(size,-1);
            for(int n=0;n<size;n++){
                if(range[n] == 1)
                {
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            //cout << "size m" << m << endl;
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0));
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else
        {
            cout << os_map[k] << endl;
            cout << label << endl;
            cout << "Contact error!!!" << endl;
            exit(1);
        }

    }
    for (int k =0 ;k <os_prf.size();k++)
    {
        istringstream is(os_prf[k]);
        string label;
        is >> label;
        if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
        else{
            cout << "PRF error!!!" << endl;
            exit(1);
        }

    }
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}

//vec_int Map::load_data (string file, int sep_cutoff, vec_bool &range){
vec_int Map::load_data (string file, string &seq, int sep_cutoff){
    vec_int n2m;
    vec_int m2n;
    string line;
    ifstream in(file);
    vec_bool range;
    while(getline(in,line)){
        istringstream is(line);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE"){
            int size; is >> size;
            
            //if(range.size() > 0){range.resize(size,0);} // if range previously defined fill rest with "0"
            //else{ range.resize(size,1);} // else set all to "1"
            range.resize(size,1);
            
            int m = 0;n2m.resize(size,-1);
            for(int n = 0; n < size; n++){
                if(range[n] == 1){
                    n2m[n] = m;
                    m2n.push_back(n);
                    m++;
                }
            }
            
            mtx.resize(m,vector<double>(m,0));
            prf.resize(m,vector<double>(4,0)); // modified by HX at 2022.02.18
            aa.resize(m,'N');
            ss.resize(m,'N');
        }
        else if(label == "CON"){
            int i, j; is >> i >> j;
            if(range[i] == 1 and range[j] == 1){
                double sco;
                if(abs(j-i) >= sep_cutoff){
                    if(is >> sco){}else{sco = 1;}
                    mtx[n2m[i]][n2m[j]] = sco;
                    mtx[n2m[j]][n2m[i]] = sco;
                }
            }
        }
        else if(label == "PRF"){
            int i; is >> i;
            if(range[i] == 1)
            {
                char tmp;
                is >> tmp; aa[n2m[i]] = tmp;
				seq += tmp;
                is >> tmp; ss[n2m[i]] = tmp;
                double val;
                int j = 0;
                while(is >> val){prf[n2m[i]][j] = val;j++;}
            }
        }
    }
    in.close();
    for(int i=0; i < mtx.size(); i++){
        vec_i.push_back(vector<int>());
        for(int j=0; j < mtx.size(); j++){
            if(i == j){
                if(vec_i[i].empty()){vec_div.push_back(0);}
                else{vec_div.push_back(vec_i[i].size());}
            }
            if(mtx[i][j] > 0){
                vec_i[i].push_back(j);
            }
        }
        if(vec_i[i].size() > 0){vec.push_back(i);}
    }
    M2N = m2n;
    return m2n;
}

int Map::printout(string outfile){
    ofstream out;
    out.open(outfile);
    for (int k=0;k<os_map.size();k++)
    {
        out << os_map[k];
    }
    for (int k =0 ;k<os_prf.size();k++){
        out << os_prf[k];
    }
    out.close();
    return 0;
}

void mod_gap(vec_double &gap, vec_char &ss, double &gap_ss_w){
    for(int i = 0; i < ss.size()-1; i++){
        //if((ss[i] == 'H' && ss[i+1] == 'H') || (ss[i] == 'E' && ss[i+1] == 'E')){gap[i] *= gap_ss_w;}
        if((ss[i] == 'P' && ss[i+1] == 'P') || (ss[i] == 'P' && ss[i+1] == 'P')){gap[i] *= gap_ss_w;}
    }
}
