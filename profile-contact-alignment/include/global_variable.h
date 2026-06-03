/*************************************************************************
    > File Name: global_variable.h
    > Author: hongxu
    > Mail: hongxu@hust.edu.cn 
    > Created Time: Sun May 15 16:19:47 2022
 ************************************************************************/

#ifndef __GLOBAL_VARIABLE_H
#define __GLOBAL_VARIABLE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include <algorithm>
#include <stack>
#include <tuple>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <map>
#include <getopt.h>

using namespace std;

typedef vector<vector<int> > mtx_int;
typedef vector<vector<double> > mtx_double;
typedef vector<int> vec_int;
typedef vector<double> vec_double;
typedef vector<char> vec_char;
typedef vector<bool> vec_bool;
typedef vector<string> vec_string;

typedef vector<int> vec_int;
typedef Eigen::Matrix<double,4,4> Matrix4d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;

//string& replace_all(string& src, const string& old_value, const string& new_value);



inline bool exists (const std::string& name) {
    ifstream f(name);
    return f.good();
}

inline float HammingDistance(string seq1,string seq2){
    int count = 0;
    if (seq1.length() != seq2.length()) {cout << "Hamming the length of seq1 and seq2 are not equal!!";exit(1);}
    for(int i=0;i<seq1.length();i++){
        if (seq1[i] != seq2[i]) count++;
    }
    return 1.0 * count / seq1.length();
}

inline string& replace_all(string& src, const string& old_value, const string& new_value){
    // replace "T" ----> "U"
    for (string::size_type pos(0); pos != string::npos; pos += new_value.length()) {
        if ((pos = src.find(old_value, pos)) != string::npos){
            src.replace(pos, old_value.length(), new_value);
        }
        else break;
    }
    return src;
}


inline bool sstest(string ss){
    // (<[{ABCD dcba}]>)
    int l1 = count(ss.begin(),ss.end(),'(');
    int l2 = count(ss.begin(),ss.end(),'<');
    int l3 = count(ss.begin(),ss.end(),'[');
    int l4 = count(ss.begin(),ss.end(),'{');
    int l5 = count(ss.begin(),ss.end(),'A');
    int l6 = count(ss.begin(),ss.end(),'B');
    int l7 = count(ss.begin(),ss.end(),'C');
    int l8 = count(ss.begin(),ss.end(),'D');
    int r1 = count(ss.begin(),ss.end(),')');
    int r2 = count(ss.begin(),ss.end(),'>');
    int r3 = count(ss.begin(),ss.end(),']');
    int r4 = count(ss.begin(),ss.end(),'}');
    int r5 = count(ss.begin(),ss.end(),'a');
    int r6 = count(ss.begin(),ss.end(),'b');
    int r7 = count(ss.begin(),ss.end(),'c');
    int r8 = count(ss.begin(),ss.end(),'d');
    if(l1 != r1 || l2 != r2 || l3 != r3 || l4 != r4 || l5 != r5 || l6 != r6 || l7 != r7 || l8 != r8)
        return false;
    else
        return true;
}

#endif