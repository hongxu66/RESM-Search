/*************************************************************************
    > File Name: load_data.cpp
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Thu Oct 20 11:06:43 2022
 ************************************************************************/

#include "load_data.h"

void Options::get_opt(vec_string opt){
    ////////////////////////////////////////////////////////////////////////////////
    auto require_value = [&](int index, const string &name) -> const string& {
        if(index + 1 >= (int)opt.size() || opt[index + 1].empty()){
            cout << "ERROR: " << name << " requires a value\n";
            exit(1);
        }
        const string &value = opt[index + 1];
        bool negative_number = value.size() > 1 && value[0] == '-' && ((value[1] >= '0' && value[1] <= '9') || value[1] == '.');
        if(value.substr(0,1) == "-" && !negative_number){
            cout << "ERROR: " << name << " requires a value\n";
            exit(1);
        }
        return value;
    };
    for (int a = 0; a < opt.size(); a++)
    {
        string arg = opt[a];
        if (arg.substr(0,1) == "-")
        {
            if(arg == "-a"){filea = require_value(a,arg); a++;}
            else if(arg == "-b"){fileb = require_value(a,arg); a++;}
//            else if(arg == "-range_a"){
//                while(a+1 < opt.size() && opt[a+1].substr(0,1) != "-")
//                {
//                    string r = opt[a+1];
//                    int i = stoi(r.substr(0,r.find('-')));
//                    int j = stoi(r.substr(r.find('-')+1));
//
//                    if(j+1 > range_a.size()){range_a.resize(j+1,0);}
//                    for(int n = i; n <= j; n++){range_a[n] = 1;}
//                    a++;
//                }
//            }
//            else if(arg == "-range_b"){
//                while(a+1 < opt.size() && opt[a+1].substr(0,1) != "-")
//                {
//                    string r = opt[a+1];
//                    int i = stoi(r.substr(0,r.find('-')));
//                    int j = stoi(r.substr(r.find('-')+1));
//
//                    if(j+1 > range_b.size()){range_b.resize(j+1,0);}
//                    for(int n = i; n <= j; n++){range_b[n] = 1;}
//                    a++;
//                }
//            }
            else if(arg == "-use_prf"){use_prf = true;}
            else if(arg == "-prf_w"){prf_w = stod(require_value(a,arg)); a++;}
            else if(arg == "-use_gap_ss"){use_gap_ss = true;}
            else if(arg == "-gap_ss_w"){gap_ss_w = stod(require_value(a,arg)); a++;}
            else if(arg == "-gap_o"){gap_open = stod(require_value(a,arg)); a++;}
            else if(arg == "-gap_e"){gap_ext  = stod(require_value(a,arg)); a++;}
            else if(arg == "-sep_cut"){sep_cutoff  = stoi(require_value(a,arg)); a++;}
            else if(arg == "-iter"){iter  = stoi(require_value(a,arg)); a++;}
            else if(arg == "-silent"){silent = true;}
            else if(arg == "-prob"){prob = stof(require_value(a,arg));a++;}
            else if(arg == "-msa") {msa = true;}
            else if(arg == "-gen") {gen = true;}
            else if(arg == "-threads" || arg == "-num_threads") {numThreads = stoi(require_value(a,arg));a++;}
            else if(arg == "-lengthC") {lengthC = stoi(require_value(a,arg));a++;}
            else if(arg == "-lengthE") {lengthE = stoi(require_value(a,arg));a++;}
            else if(arg == "-overlap") {overlap = stoi(require_value(a,arg));a++;}
            else if(arg == "-msafile") {msafile = require_value(a,arg);a++;}
            else if(arg == "-dbnfile") {dbnfile = require_value(a,arg);a++;}
            else if(arg == "-stofile") {stofile = require_value(a,arg);a++;}
            else if(arg == "-msafileb") {msafileb = require_value(a,arg);a++;}

        }
    }

    if(exists(filea) == 0 || exists(fileb) == 0 || filea.empty() || fileb.empty())
    {
        cout << "-------------------------------------------------------------------\n";
        cout << "                          CMAP_SEARCH                              \n";
        cout << "-------------------------------------------------------------------\n";
        cout << "  -a             RNA secondary structure file            [REQUIRED]\n";
        cout << "  -b             RNA sequence file                       [REQUIRED]\n";
        cout << "  -gap_o         gap opening penalty         [Default=" << gap_open << "]\n";
        cout << "  -gap_e         gap extension penalty       [Default=" << gap_ext << "]\n";
        cout << "  -sep_cut       seq seperation cutoff       [Default=" << sep_cutoff << "]\n";
        cout << "  -iter          number of iterations        [Default=" << iter << "]\n";
        cout << "  -threads       OpenMP worker threads       [Default=" << numThreads << "]\n";
        cout << "  -silent        \n";
        cout << "-------------------------------------------------------------------\n";
        cout << " Advanced options\n";
        cout << "-------------------------------------------------------------------\n";
        cout << "  -range_a       trim map A to specified range(s) (eg. 0-20 50-100)\n";
        cout << "  -range_b       trim map B to specified range(s)\n";
        cout << "-------------------------------------------------------------------\n";
        cout << " Experimental features\n";
        cout << "-------------------------------------------------------------------\n";
        cout << "  -use_gap_ss    penalize gaps at secondary structure elements(SSE)\n";
        cout << "  -gap_ss_w      gap penality weight at SSE  [Default=" << gap_ss_w << "]\n";
        cout << "  -use_prf       use sequence profile\n";
        cout << "  -prf_w         profile weight              [Default=" << prf_w << "]\n";
        cout << "  -prob          sequence profile initialization       [Default=" << prob << "]\n";
        cout << "  -mul           single dbn mapped to multiple sequences   [Default=" << msa << "]\n";
        cout << "-------------------------------------------------------------------\n";
//        exit(1);
    }
    else if(silent == false)
    {
        cout << "OPT -------------------------------------------------------------------\n";
        cout << "OPT                          RESM_Search                               \n";
        cout << "OPT -------------------------------------------------------------------\n";
        cout << "OPT   -a          " << filea       << endl;
        cout << "OPT   -b          " << fileb       << endl;
        cout << "OPT   -gap_o      " << gap_open     << endl;
        cout << "OPT   -gap_e      " << gap_ext      << endl;
        cout << "OPT   -sep_cut    " << sep_cutoff   << endl;
        cout << "OPT   -iter       " << iter         << endl;
        cout << "OPT   -threads    " << numThreads   << endl;
        cout << "OPT   -silent     " << silent       << endl;
        cout << "OPT -------------------------------------------------------------------\n";
        cout << "OPT   -use_gap_ss  " << use_gap_ss  << endl; if(use_gap_ss == true){cout << "OPT   -gap_ss_w    " << gap_ss_w << endl;}
        cout << "OPT   -use_prf     " << use_prf     << endl; if(use_prf    == true){cout << "OPT   -prf_w       " << prf_w    << endl;}
        cout << "OPT -------------------------------------------------------------------\n";
    }
    if(gap_open > 0 || gap_ext > 0)
    {
        cout << "ERROR: gap penality should be < 0\n";
        exit(1);
    }
    if(numThreads < 1)
    {
        cout << "ERROR: -threads should be >= 1\n";
        exit(1);
    }
    if(lengthC < 1 || lengthE < 1 || overlap < 0 || overlap >= lengthE)
    {
        cout << "ERROR: require lengthC >= 1, lengthE >= 1 and 0 <= overlap < lengthE\n";
        exit(1);
    }
}

void readdbnfile(vec_string &seqnames, vec_string &seqs, vec_string &seqss,const std::string filename){
    ifstream in(filename.c_str());
    //seqdbn *queryseq = new seqdbn;
    if(!in)
    {
        std::cout << "dbnfile is empty!!!" << std::endl;
        exit(1);
    }
    else{
        string s;
        string temseq,temid,temss;
        int i=0;
        while (getline(in,s)){
			if(s.empty()) continue;
            i++;
            cout << s << endl;
            //std::cout << "fasta file is " << s << std::endl;
            if(i%3==1)
            {
                temid = s;
                seqnames.push_back(temid.substr(1));
            }
            if(i%3==2)
            {
                string s1 = replace_all(s,"T","U");
                temseq = s1;
                seqs.push_back(temseq);
            }
            if(i%3==0)
            {
                temss = s;
                seqss.push_back(temss);
            }
        }
    }
    cout << "read " << endl;
    cout << seqnames.size() << endl;
    cout << seqs.size() << endl;
    cout << seqss.size() << endl;
    for (int i = 0;i<seqnames.size();i++){
        cout << i << " " << seqnames[i] << " " << seqss[i] << " " << seqs[i] << endl;
        if (seqss[i].length() != seqs[i].length()){
            cout << "DBN fasta is not equal to ss " << endl;
            exit(1);
        }
    }
    cout << "read dbn done" << endl;
    in.close();
}

void readfastafile(vec_string & refnames, vec_string & refseqs, const std::string filename){
    //  read fasta file inincluding multiple sequences
    ifstream in(filename.c_str());
    if(!in)
    {
        std::cout << "fastafile is empty!!!" << std::endl;
        exit(1);
    }
    else{
        std::string line, name, seq;
        while( std::getline( in, line ) ){
			if (line.empty()) continue;
            if( line.empty() || line[0] == '>' ){
                if( !name.empty() ){
                    //std::cout << name << " : " << seq << std::endl;
                    refnames.push_back(name);
                    refseqs.push_back(replace_all(seq,"T","U"));
                    name.clear();
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
            refnames.push_back(name);
            refseqs.push_back(replace_all(seq,"T","U"));
//			cout << seq << endl;
        }
    }
    if (refnames.size() != refseqs.size())
    {
        cout << "The size of Queryid is not equal to the size of Queryseqs\n";
        exit(1);
    }
    cout << "read fasta done " << endl;
    in.close();
}

void readfastafile(Options OPT, vec_string & refnames, vec_string & refseqs, const std::string filename){
    //  read fasta file inincluding multiple sequences
	cout << "filename:" << filename << endl;
    ifstream in(filename.c_str());
    if(!in)
    {
        std::cout << "fastafile is empty!!!" << std::endl;
        exit(1);
    }
    else{
        std::string line, name, seq;
        while( std::getline( in, line ) ){
			if (line.empty()) continue;
            if( line.empty() || line[0] == '>' ){
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
                            refnames.push_back(newname);
                            refseqs.push_back(replace_all(temseq,"T","U"));
                        }
                    }
                    else{
                        refnames.push_back(name);
                        refseqs.push_back(replace_all(seq,"T","U"));
//                        name.clear();
                    }
                    name.clear();
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
            if(seq.length() >= OPT.lengthC){
                        int seqLength = seq.length();
                        int K = ceil(1.0* (seqLength - OPT.lengthE) / (OPT.lengthE - OPT.overlap)) + 1;
                        for(int ii=0;ii<K;ii++){
                            string temseq;
                            string newname = name + ":" + to_string(ii);
                            if(ii == K -1 ) temseq = seq.substr((OPT.lengthE - OPT.overlap)  * ii);
                            else
                                temseq = seq.substr((OPT.lengthE - OPT.overlap) * ii,OPT.lengthE);
                            refnames.push_back(newname);
                            refseqs.push_back(replace_all(temseq,"T","U"));
                        }
                    }
                    else{
                        refnames.push_back(name);
                        refseqs.push_back(replace_all(seq,"T","U"));
//                        name.clear();
                    }
                    name.clear();
        }
    }
    if (refnames.size() != refseqs.size())
    {
        cout << "The size of Queryid is not equal to the size of Queryseqs\n";
        exit(1);
    }
    in.close();
    cout << "read fasta done " << endl;
}

//sequence * readsinglefasta(const std::string filename);

seqdbn * readdbnfile(const std::string filename){
    // read fasta file then return id and sequence
    ifstream in(filename.c_str());
    seqdbn *queryseq = new seqdbn;
    if(!in)
    {
        std::cout << "fastafile is empty!!!" << std::endl;
        exit(1);
    }
    else{
        string s;
        int i=0;
        while (getline(in,s)){
            i++;
            //std::cout << "fasta file is " << s << std::endl;
            if(i==1)
                queryseq->id = s;
            if(i==2)
            {
                string s1 = replace_all(s,"T","U");
                queryseq->seq = s1;
            }
            if(i==3)
            {
                queryseq->ss = s;
                break;
            }
        }
    }
    //cout << queryseq->id << " " << queryseq->seq << endl;
    in.close();
    return queryseq;

}

string readctfile(const std::string filename){
    // ct file
    // read dbn file, dbn file includes sequence and secondary structure
    ifstream in(filename.c_str());
    string ss;
    if(!in)
    {
        std::cout << "dbnfile file is empty!!!" << std::endl;
        exit(1);
    }
    else{
        string s;
        int i=0;
        while (getline(in,s)){
            i ++;
            std::cout << "dbnfile is " << s << std::endl;
            if (i==3)
                ss = s;
        }
    }
    in.close();
    return ss;

}


sequence * readsinglefasta(const std::string filename){
    // read fasta file then return id and sequence
    ifstream in(filename.c_str());
    sequence *queryseq = new sequence;
    string seq;
    if(!in)
    {
        std::cout << "fastafile is empty!!!" << std::endl;
        exit(1);
    }
    else{
        string s;
        int i=0;
        while (getline(in,s)){
            i++;
            //std::cout << "fasta file is " << s << std::endl;
            if(i==1)
                queryseq->id = s;
            else
            {
                seq += s;
                //break;
            }
        }
    }
    queryseq->seq = replace_all(seq,"T","U");
    //cout << queryseq->id << " " << queryseq->seq << endl;
    in.close();
    return queryseq;

}
