/*
 * Kinpute is a program to do imputation of sequence genotype data based on
 * identity by descent information.
 * 
 * Copyright 2018 Aisha El Sherbiny and Mark Abney
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */ 

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <functional>
#include <sstream>
#include <cmath>
#include <time.h>
#include <omp.h>
#include <string.h>
//#include "mapsampler_api.h"
using namespace std;
struct Key
{
  string r1;
  string r2;
  int snp;
  bool operator==(const Key &other)const
  {return (r1 == other.r1
           && r2 == other.r2
           && snp == other.snp);
  }
};
namespace std {
template<>
struct hash<Key>
{
  std::size_t operator()(const Key&k)const
  {
    using std::size_t;
    using std:: hash;
    using std::string;
    return((hash<string>()(k.r1)
            ^(hash<string>()(k.r2) << 1)) >>1)
        ^(hash<int>()(k.snp)<< 1);
  }
};
        
}
typedef pair<string,string> ID_pair;
struct hash_name {
  size_t operator()(const ID_pair &name ) const
  {
    return hash<string>()(name.first) ^ hash<string>()(name.second);
  }
};
typedef pair<string,int> US_pair;
struct hash_name1 {
  size_t operator()(const US_pair &name ) const
  {
    return hash<string>()(name.first) ^ hash<int>()(name.second);
  }
};
void Upr2(vector<double>& v,unordered_map <US_pair, vector<double>, hash_name1>& Upr, double P, double gen, double fu, int snp, string U);
void get_delta(vector<double>& delta1, vector<double>& delta2, vector<double>& delta3, vector<double>& delta4, vector<double>& delta5, vector<double>& delta6,
               vector<double>& delta7, vector<double>& delta8, vector<double>& delta9,
               vector<double>& v, int Ns1);

void get_W(vector<double>& W0, vector<double>& W1, vector<double>& W2, double delta1,
           double delta2, double delta3, double delta4, double delta5, double delta6,
           double delta7, double delta8, double delta9, double p);
void prGUS(int S, double gen, unordered_map<US_pair, vector<double>, hash_name1>& Upr, int snp, string U, double prG1G2S12[][3][3]);
void vec_mat_mult(double e[],double A[][9],int m, int n,double eA[]);
void max_sec(vector<double>& v3, vector<double>& max, int N);

int main (int argc, char* argv[]) {
  string myfile1;
  string myfile2;
  string myfile4;
  string myfile5;
  string myfile6;
  string myfile7;
  string myfile8;
  bool flag = false;
  if(argc < 13){
    cout<<"Usage is: ./" << argv[0] << " -ibd <infile1> -map <infile2> -u <infile3> -r <infile4> -seq <infile5> -prior <infile6> -o <outfile>" << endl;
    exit(0);
  }else{
    for(int i=1; i<argc; i++){
      if(i+1 != argc){
        if(strcmp(argv[i], "-ibd")==0){
          myfile1 = argv[i+1];
        }else if(strcmp(argv[i], "-map")==0){
          myfile2 = argv[i+1];
        }else if(strcmp(argv[i], "-u")==0){
          myfile4 = argv[i+1];
        }else if(strcmp(argv[i], "-r")==0){
          myfile5 = argv[i+1];
        }else if(strcmp(argv[i], "-seq")==0){
          myfile6 = argv[i+1];
        }else if(strcmp(argv[i], "-prior")==0){
          myfile7 = argv[i+1];
        }else if(strcmp(argv[i], "-o")==0){
          myfile8 = argv[i+1];
          // }else if(strcmp(argv[i], "-f")==0){
          //   istringstream(argv[i+1]) >> flag;
        }
      }
    }
  }
  string st1;
  ifstream ifile04(myfile5);
  if(! ifile04){
    cerr << "Error: Could not find file " << myfile5 << "\n";
  } 
  vector<string> SID;                              
  while(ifile04 >> st1){
    SID.push_back(st1);
  }
  ifstream ifile05(myfile4);
  if(! ifile05){
    cerr << "Error: Could not find file " << myfile4 << "\n";
  } 
  vector<string> UID;                              
  while(ifile05 >> st1){
    UID.push_back(st1);
  }

  int val1, val2;
  string str1, str2, str3, str4;
  string line0;
  int val3, val4;
  vector<int>comm;
  // vector<int> pos90;
  // ifstream infile00(myfile6);
  // if(! infile00){
  //   cerr << "Error: Could not find file " << myfile6 << "\n";
  // } 
  // while(getline(infile00,line0)){
  //   stringstream s000(line0);
  //   s000 >> val1 >> val2 >> val3 >> val4;
  //   vector <string> v1;
  //   string s;
  //   while(s000 >> s){
  //     v1.push_back(s);
  //   }
  //   if(find(pos90.begin(), pos90.end(), val4)==pos90.end()){
  //     pos90.push_back(val4);
  //   }
  // }

  string line1;
  int pos ;
  double fre1, fre2, fre3, fre4, fre5, fre6;
  string snpid, allele00, allele1, bi;
  unordered_map <int, vector<string> > allele ;
  unordered_map <int, vector<string> > allele0 ;
  unordered_map<US_pair, vector<double>, hash_name1> gpro1;
  vector<int> pos2;

  ifstream infile070(myfile6);
  if(! infile070){
    cerr << "Error: Could not find file " << myfile6 << "\n";
  }

  if (!myfile7.empty())
    flag = true;
  else
    cerr << "Warning: No genotype prior probabilities. Using allele frequencies instead.\n";
  
  if(flag) {
    while (getline(infile070,line1)){
      stringstream s8(line1);
      s8 >> val1 >> val2 >> val3 >> val4;
      vector <string> v1;
      string s;
      while(s8 >> s){
        v1.push_back(s);
      }
      if(find(pos2.begin(), pos2.end(), val4)==pos2.end()){
        pos2.push_back(val4);
      }

      for(int k = 0; k < v1.size(); k += 2){
        str1 = v1[k];
        str2 = v1[k+1];
        if(str1 != str2){
          break;
        }
      }
      if(str1 == str2) {
        for(int k = 0; k < v1.size(); k++){
          str1 = v1[k];
          if(str1 != "0"){
            break;
          }
        }
        for(int k = 0; k < v1.size(); k++) {
          str2 = v1[k];
          if(str2 != str1 && str2 != "0") {
            break;
          } else { // MA: Prevent str2=0 if last genotype is 0 0 and only one allele exists
            str2 = str1;
          }
        }
      }
      vector<string> v2;
      v2.push_back(str1);
      v2.push_back(str2);
      allele0[val4] = v2;
    }
  }

  ifstream infile07(myfile6);
  if(! infile07){
    cerr << "Error: Could not find file " << myfile6 << "\n";
  } 
  if(!flag) {
    while (getline(infile07,line1)){
      stringstream s8(line1);
      s8 >> val1 >> val2 >> val3 >> val4;
      vector <string> v1;
      string s;
      while(s8 >> s){
        v1.push_back(s);
      }
      if(find(pos2.begin(), pos2.end(), val4)==pos2.end()) {
        pos2.push_back(val4);
        comm.push_back(val4);
        for(int k = 0; k < v1.size(); k += 2) {
          str1 = v1[k];
          str2 = v1[k+1];
          if(str1 != str2){
            break;
          }
        }
        if(str1 == str2) {
          for(int k = 0; k < v1.size(); k++) {
            str1 = v1[k];
            if(str1 != "0"){
              break;
            }
          }
          for(int k = 0; k < v1.size(); k++) {
            str2 = v1[k];
            if(str2 != str1 && str2 != "0") {
              break;
            } else { // MA: Prevent str2=0 if last genotype is 0 0 and only one allele exists
              str2 = str1;
            }
          }
        }
        vector<string> v2;
        v2.push_back(str1);
        v2.push_back(str2);
        allele[val4] = v2;
        string str3;
        string str4;
        for(int i=0; i<v1.size(); i+=2){
          str3 = v1[i];
          str4 = v1[i+1];
          vector<double>v;
          if(str3 == str4 && str3 == str1){
            v.push_back(1);
            v.push_back(0);
            v.push_back(0);
          }else if(str3 ==str4 && str3 == str2){
            v.push_back(0);
            v.push_back(0);
            v.push_back(1);
          }else if(str3 !=str4){
            v.push_back(0);
            v.push_back(1);
            v.push_back(0);
          }else if(str3 =="0" && str4 =="0"){
            v.push_back(0);
            v.push_back(0);
            v.push_back(0);
          }
          gpro1.insert(make_pair(make_pair(SID[i/2],val4),v));
        }
      }
    }
  }

  string line;
  string str,str5;
  int val;
  double val11, val22, val33, val44, val55, val66, val77, val88, val99, val100;
  ifstream ifile09(myfile1);
  if(! ifile09){
    cerr << "Error: Could not find file " << myfile1 << "\n";
  } 
  getline(ifile09,line1);
  stringstream ss00(line1);
  string strr0;
  vector<string> vv0;
  while(ss00 >> strr0){
    vv0.push_back(strr0);}
  int Ns=stoi(vv0[3]);
  string chro = vv0[2];
  // Require that the IBD file have all nine identity coefficients for each SNP.
  if (stoi(vv0[0]) != 9) {
    cerr << "The IBD file does not have the correct number of identity coefficients.\n";
    exit(2);
  }
  
  vector<int> pos9;
  vector<int> priv;
  map<pair<string, int>, vector<double>> PGu;
  if(flag) {
    ofstream output22("ppf");
    ifstream infile10(myfile7);
    if(! infile10){
      cerr << "Error: Could not find file " << myfile7 << "\n";
    } 
    string line2;
    string  str6, str7, str8;
    getline(infile10,line1);
    stringstream s1(line1);
    s1 >> str1 >> str2 >> val1 >> str3 >> str4;

    while(getline(infile10,line2)) {
      bool f = true;
      stringstream s2(line2);
      s2 >> str5 >> str6 >> val2 >> str7 >> str8;
      if(val2 > val1) {
        istringstream s3(line1);
        s3 >> str1 >> str2 >> val1 >> str3 >> str4;
        vector<double> v1;
        double val22;
        while(s3 >> val22){
          v1.push_back(val22);
        }
        if(find(pos2.begin(), pos2.end(), val1) != pos2.end()) {
          string al0 = allele0[val1][0];
          string al1 = allele0[val1][1];
          // if(al0 != al1) {
          // if((str3 == al0 || str3 == al1) && (str4 == al0 || str4 == al1)) {
          if((str3 == al0 || str4 == al0) && (str3 == al1 || str4 == al1)) {
            vector<string> v3;
            v3.push_back(str3);
            v3.push_back(str4);
            allele[val1] = v3;
            pos9.push_back(val1);
            output22<<"---"<<" "<<chro<<":"<<val1<<":"<<v3[0]<< 
                ":"<<v3[1]<<" "<<val1<<" "<<v3[0]<<" "<<
                v3[1]<<" ";
            for(int i=0; i<v1.size(); i++)
              output22 << v1[i] << " ";
            output22 << endl;
          }
          // }
        }
      } else {
        f=false;
        getline(infile10,line1);
        stringstream s1(line1);
        s1 >> str1 >> str2 >> val1 >>str3 >>str4;
      }
      if(f){
        line1 = line2;
        val1 = val2;
      }
    }
    // Still need to process the last line:
    istringstream s3(line1);
    s3 >> str1 >> str2 >> val1 >> str3 >> str4;
    vector<double> v1;
    double val22;
    while(s3 >> val22){
      v1.push_back(val22);
    }
    if(find(pos2.begin(), pos2.end(), val1) != pos2.end()) {
      string al0 = allele0[val1][0];
      string al1 = allele0[val1][1];
      // if(al0 != al1) {
      // if((str3 == al0 || str3 == al1) && (str4 == al0 || str4 == al1)) {
      if((str3 == al0 || str4 == al0) && (str3 == al1 || str4 == al1)) {
        vector<string> v3;
        v3.push_back(str3);
        v3.push_back(str4);
        allele[val1] = v3;
        pos9.push_back(val1);
        output22<<"---"<<" "<<chro<<":"<<val1<<":"<<v3[0]<< 
            ":"<<v3[1]<<" "<<val1<<" "<<v3[0]<<" "<<
            v3[1]<<" ";
        for(int i=0; i<v1.size(); i++)
          output22 << v1[i] << " ";
        output22 << endl;
      }
      // }
    }
  }

  if(!flag)
    pos9 = pos2;
  ifstream infile077(myfile6);
  if(! infile077){
    cerr << "Error: Could not find file " << myfile6 << "\n";
  }
  
  if(flag) {
    while (getline(infile077,line1)) {
      stringstream s8(line1);
      s8 >> val1 >> val2 >> val3 >> val4;
      vector <string> v1;
      string s;
      while(s8 >> s) {
        v1.push_back(s);
      }
      string str5, str6;
      for(int k = 0; k < v1.size(); k += 2) {
        str5 = v1[k];
        str6 = v1[k+1];
        if(str5 != str6){
          break;
        }
      }
      if(str5 == str6) {
        for(int k = 0; k < v1.size(); k++) {
          str5 = v1[k];
          if(str5 != "0"){
            break;
          }
        }
        for(int k = 0; k < v1.size(); k++) {
          str6 = v1[k];
          if(str6 != str5 && str6 != "0"){
            break;
          } else { // MA: Prevent str6=0 if last genotype is 0 0 and only one allele exists
            str6 = str5;
          }
        }
      }

      if(find(pos9.begin(), pos9.end(), val4) != pos9.end()) {
        string str1 = allele[val4][0];
        string str2 = allele[val4][1];
        // if((str1 == str5 || str1 == str6) && (str2 == str5 || str2 == str6)) {
        if((str5 == str1 || str5 == str2) && (str6 == str1 || str6 == str2)) {
          string str3;
          string str4;
          for(int i = 0; i < v1.size(); i += 2){
            str3 = v1[i];
            str4 = v1[i+1];
            vector<double> v;
            if(str3 == str4 && str3 == str1){
              v.push_back(1);
              v.push_back(0);
              v.push_back(0);
            }else if(str3 == str4 && str3 == str2){
              v.push_back(0);
              v.push_back(0);
              v.push_back(1);
            }else if(str3 != str4){
              v.push_back(0);
              v.push_back(1);
              v.push_back(0);
            }else if(str3 == "0" && str4 =="0"){
              v.push_back(0);
              v.push_back(0);
              v.push_back(0);
            }
            gpro1.insert(make_pair(make_pair(SID[i/2],val4),v));
          }
        }
      }
    }
  }          



  unordered_map <int, vector<double> > freq1;
  unordered_map <int, double > maf;
  unordered_map <int, int > mac1;
  vector<int>mac;
  for(int& pos:pos9){
    vector<string>v2;
    vector<string>v;
    v2 = allele[pos];
    string str1 = v2[0];
    string str2 = v2[1];
    for(string& str:SID){
      vector<double>v1;
      v1 = gpro1[make_pair(str,pos)];
      if(v1.size() ==0){
        cout<<v1.size()<<" "<<str<<" "<<pos<<endl;
        exit(1);
      }
      auto l1 = max_element(begin(v1), std::end(v1));
      double d2 = distance(begin(v1),l1);
      if(d2==0){
        v.push_back(str1);
        v.push_back(str1);
      }else if(d2==1){
        v.push_back(str1);
        v.push_back(str2);
      }else if(d2==2){
        v.push_back(str2);
        v.push_back(str2);
      }
    }
    vector<double>v3;
    double count0 = 1;
    double count1 = 1;
    for(int i=0; i<v.size(); i++){
      int found =0;
      for(int j=0; j<i; j++){
        if(v[i] == v[j]) found++;
      }
      if(found ==0){
        for(int j=i+1; j<v.size(); j++){
          if(v[i] == v[j] && v[j]==str1 && v[i] !="0" && v[j]!="0")
          {count0++;}else if(v[i]==v[j] && v[j]== str2 && v[i] !="0" && v[j]!="0"){
            count1++;}
        }
      }
    }
    if(str1 ==str2 && count0>0 ){
      count1=0;
    }else if(str1 ==str2 && count1>0 ){
      count0=0;
    }

    double countot = count0 + count1;
    v3.push_back(count0/countot);
    v3.push_back(count1/countot);
    freq1[pos]=v3;
    if((count0/countot)<= .5){
      maf[pos]= count0/countot;
    }else{
      maf[pos]= count1/countot;
    }
  }
  ifstream ifile07(myfile2);
  if(! ifile07){
    cerr << "Error: Could not find file " << myfile2 << "\n";
  } 
  vector<string> snpID;
  vector<int> pos3;
  while(ifile07 >> val11 >>st1 >>val22 >>val3){
    snpID.push_back(st1);
    pos3.push_back(val3);
  }


  ifstream ifile08(myfile1);
  if(! ifile08){
    cerr << "Error: Could not find file " << myfile1 << "\n";
  } 
  getline(ifile08,line1);
  stringstream sss(line1);
  sss >> str1 >> str2 >> str3 >> str4;
  string st0;
  vector<string> v0;
  while(sss >> st0){
    v0.push_back(st0);
  }
  vector<int>posibd0;
  for(int i=0; i<snpID.size(); i++){
    if(find(v0.begin(), v0.end(), snpID[i]) !=v0.end()){
      posibd0.push_back(pos3[i]);
    }
  }


  unordered_map <int, bool> genbol1 ;
  unordered_map <int, int> interval ;
  unordered_map <int, int> posmap ;

  vector<int>posibd;
  int kk=0;
  for(int i=0; i<posibd0.size(); i++){
    if(find(pos9.begin(), pos9.end(), posibd0[i]) !=pos9.end()){
      posibd.push_back(posibd0[i]);
      posmap[i] = kk;
      kk++;
      genbol1.insert(make_pair(i,true));
    }else{
      genbol1.insert(make_pair(i,false));
    }
  }
  double countpos9=0.0;
  double cpos91=0.0;
  double cpos92=0.0;
  double cpos93=0.0;
  unordered_map <int, vector<int>> interval1 ;
  for(int i=0; i<pos9.size(); i++){
    if(posibd0[0]>pos9[i]){
      interval1[0].push_back(i);
      cpos91++;
      countpos9++;
    }
  }
  for(int j=0; j<posibd0.size(); j++){
    vector<int>v;
    if(j>posibd0.size()-2) break;
    for(int i=0; i<pos9.size(); i++){
      if(pos9[i] >=posibd0[j] && pos9[i]<posibd0[j+1]){
        if(j==0){
          interval1[0].push_back(i);
          countpos9++;
          cpos92++;
        }else{
          v.push_back(i);
        }
      }
    }
    if(j>0){
      interval1.insert(make_pair(j,v));
      countpos9= countpos9+ v.size();
      cpos92= cpos92+ v.size();
    }
  }
  for(int i=0; i<pos9.size(); i++){
    if(posibd0[posibd0.size()-1]<=pos9[i]){
      interval1[posibd0.size()-1].push_back(i);
      countpos9++;
      cpos93++;
    }
  }
  for(int i=0; i<pos9.size(); i++){
    for(int j=0; j<posibd0.size(); j++){
      if(j>posibd0.size()-2) break;
      if(pos9[i] >=posibd0[j] && pos9[i]<posibd0[j+1]){
        interval.insert(make_pair(i,9*j));
      }
    }
  }
  int counter1=0;
  int counter2=0;

  for(int i=0; i<pos9.size(); i++){
    if(posibd0[0]>pos9[i]){
      counter1++;}
  }
  for(int i=0; i<pos9.size(); i++){
    if(posibd0[posibd0.size()-1]<=pos9[i]){
      counter2++;}
  }
  map<pair<string, int>, double> Uf;
  unordered_map<ID_pair,vector <double>,hash_name> data;
  unordered_map <string, vector<vector <double> >> UR ;
  unordered_map <string, vector<string> > RID ;
  unordered_map <int, vector<double> > config ;
  map<pair<string, int>, int>FU;
  std::unordered_map<Key,vector<double>> data2;
  getline(ifile09,line1);
  stringstream ss0(line1);
  ss0 >> str1 >> str2 >> str3 >> str4;
  vector <double> line3;
  double value1;
  while(ss0 >>value1)
    line3.push_back(value1);
  int Nfiles = ceil(double(line3.size())/9000.0);
  double epsilon = 0.005;
  vector<string>files;
  for(int i=1; i<=Nfiles; i++){
    str1= to_string(i);
    string str = string("input")+"."+str1;
    files.push_back(str);
  }
  ifstream ifile090(myfile1);
  if(! ifile090){
    cerr << "Error: Could not find file " << myfile1 << "\n";
  } 
  getline(ifile090,line1);

  ofstream output[Nfiles];
  int l=0;
  for(int i=0; i<Ns*9; i=i+9000){
    output[l].open(files[l]);
    l++;
  }
  if(Ns <=1000){
    while(getline(ifile090,line1)){
      string str1, str2, str3, str4;
      stringstream ss(line1);
      ss >> str1 >> str2 >> str3 >> str4;
      vector <double> line;
      double value1;
      while(ss >>value1)
      {line.push_back(value1);}
      output[0]<<str1<<" "<<str2<<" "<<str3<<" "<<str4<<" ";
      for(int i=0; i<Ns*9; i++){
        output[0]<<line[i]<<" ";
      }
      output[0]<<endl;
    }
  }else{
    int m;
    while(getline(ifile090,line1)){
      string str1, str2, str3, str4;
      stringstream ss(line1);
      ss >> str1 >> str2 >> str3 >> str4;
      vector <double> line;
      double value1;
      while(ss >>value1)
        line.push_back(value1);
      int l=0;
      for(int i=0; i<line.size(); i=i+9000){
        if(i+9000 >line.size()){}
        else{
          if(i==0){
            i=i+9;
          }
          output[l]<<str1<<" "<<str2<<" "<<str3<<" "<<str4<<" ";
          for(int j=i-9; j<i+9000; j++){
            output[l]<<line[j]<<" ";
            m =j;
          }
          output[l]<<endl;
          l++;
        }
      }
      output[l]<<str1<<" "<<str2<<" "<<str3<<" "<<str4<<" ";
      for(int k=m-8; k<line.size(); k++){
        output[l]<<line[k]<<" ";
      }
      output[l]<<endl;
    }
  }
  
  vector<string>SS;
  double SScount=0; 
  ifstream ifile010(myfile1);
  if(! ifile010){
    cerr << "Error: Could not find file " << myfile1 << "\n";
  } 
  getline(ifile010,line1);
  while(getline(ifile010,line1)){
    string str1, str2, str3, str4, str5, str6;
    stringstream ss(line1);
    ss >> str1 >> str2 >> str3 >> str4;
    str5 = str1+","+str2;
    str6 = str3+","+str4;
    vector <double> line;
    double value1;
    while(ss >>value1)
    {line.push_back(value1);}
    if(find(UID.begin(), UID.end(), str5) !=UID.end()
       &&str5 == str6){
      vector<double>line2;
      for (int k=0; k<line.size(); k+=9){
        line2.push_back(line[k]);
        line2.push_back(line[k+1]);
        line2.push_back(line[k+2]);
        line2.push_back(line[k+3]);
        line2.push_back(line[k+4]);
        line2.push_back(line[k+5]);
        line2.push_back(line[k+6]);
        line2.push_back(line[k+7]);
        line2.push_back(line[k+8]);
      }
      double del1;
      for(int j=0; j<pos9.size(); j++){
        if(posibd0[0] >pos9[j]){
          del1 = line2[0];
        }else if(posibd0[posibd0.size()-1] <=pos9[j]){
          del1 = line2[line2.size()-9];
        }else{
          del1 = line2[interval[j]]+(pos9[j]-posibd0[interval[j]/9])*
              (line2[interval[j]+9]-line2[interval[j]])/(posibd0[(interval[j]/9)+1]-
                                                         posibd0[interval[j]/9]);
        }
        Uf.insert(make_pair(make_pair(str5,j),del1));
      }
    }
    if(find(SID.begin(), SID.end(), str5) !=SID.end() &&
       find(SID.begin(), SID.end(), str6) !=SID.end())
    {
      string str7 = str5+str6;
      string str8 = str6+str5;
      if(find(SS.begin(), SS.end(), str7) ==SS.end() &&
         find(SS.begin(), SS.end(), str8) ==SS.end()){
        SS.push_back(str7);
        SScount++;
      }
      int kk=0;
      for (int k=0; k<line.size(); k+=9){
        vector<double>line2;
        line2.push_back(line[k]);
        line2.push_back(line[k+1]);
        line2.push_back(line[k+2]);
        line2.push_back(line[k+3]);
        line2.push_back(line[k+4]);
        line2.push_back(line[k+5]);
        line2.push_back(line[k+6]);
        line2.push_back(line[k+7]);
        line2.push_back(line[k+8]);

        Key key1={str5,str6,kk};
        data2.emplace(key1,line2);
        Key key2={str6,str5,kk};
        data2.emplace(key2,line2);
        kk++;
      }
    }
  }
          
  if(SScount < (SID.size()*(SID.size()+1))/2){
    cout<<"IBDLD file is missing pairs fromm the reference panel. Now exiting." << endl;
    exit(1);
  }
 
  int len1 =0;
  int len2 =0;
  unordered_map<US_pair, vector<double>, hash_name1>Upr;
  unordered_map <string, bool > Rbool; 
  unordered_map <string, vector<double> > Rdelta ;
  unordered_map <string, vector<double> > URS;
  ifstream in;
  int mm=-1;
  for(int i=0; i<l; i++){
    mm++;
    UR.clear();
    RID.clear();
    in.open(files[i].c_str(), fstream::in);
    while(getline(in,line1)){
      string str1, str2, str3, str4, str5, str6;
      stringstream ss(line1);
      ss >> str1 >> str2 >> str3 >> str4;
      str5 = str1+","+str2;
      str6 = str3+","+str4;
      vector <double> line;
      double value1;
      while(ss >>value1)
        line.push_back(value1);
      if(std::find(SID.begin(), SID.end(), str5) !=SID.end()&&
         (std::find(UID.begin(), UID.end(), str6) !=UID.end())){
        vector<double>line2;
        for (int k=0; k<line.size(); k+=9){
          line2.push_back(line[k]);
          line2.push_back(line[k+1]);
          line2.push_back(line[k+2]);
          line2.push_back(line[k+3]);
          line2.push_back(line[k+4]);
          line2.push_back(line[k+5]);
          line2.push_back(line[k+6]);
          line2.push_back(line[k+7]);
          line2.push_back(line[k+8]);
        }
        len2 = line2.size()-9;
        if(mm==(l-1))
          len2 = line2.size();
        unordered_map <string, vector <vector<double> >>:: iterator itR= UR.find(str6) ;
        vector<vector <double> > v1;
        vector<string> v2;
        if(itR !=UR.end()){
          v1 = UR[str6];
          v1.push_back(line2);
          UR[str6] = v1;
          v2 = RID[str6];
          v2.push_back(str5);
          RID[str6] =v2;
        }else{
          v1.push_back(line2);
          UR.insert(make_pair(str6,v1));
          v2.push_back(str5);
          RID.insert(make_pair(str6,v2));
        }
      }else if(std::find(UID.begin(), UID.end(), str5) !=UID.end()&&
               (std::find(SID.begin(), SID.end(), str6) !=SID.end())){
        vector<double>line2;
        for (int k=0; k<line.size(); k+=9){
          line2.push_back(line[k]);
          line2.push_back(line[k+1]);
          line2.push_back(line[k+2]);
          line2.push_back(line[k+3]);
          line2.push_back(line[k+4]);
          line2.push_back(line[k+5]);
          line2.push_back(line[k+6]);
          line2.push_back(line[k+7]);
          line2.push_back(line[k+8]);
        }
        len2 = line2.size()-9;
        if(mm==(l-1))
          len2 = line2.size();
        for(int i=0; i<line2.size(); i+=9){
          double temp1 = line2[i+2];
          line2[i+2] = line2[i+4];
          line2[i+4] = temp1;
          double temp2 = line2[i+3];
          line2[i+3] = line2[i+5];
          line2[i+5] = temp2;
        }
        unordered_map <string, vector <vector<double> >>:: iterator itR= UR.find(str5) ;
        vector<vector <double> > v1;
        vector<string> v2;
        if(itR !=UR.end()){
          v1 = UR[str5];
          v1.push_back(line2);
          UR[str5] = v1;
          v2 = RID[str5];
          v2.push_back(str6);
          RID[str5] = v2;
        }else{
          v1.push_back(line2);
          UR.insert(make_pair(str5,v1));
          v2.push_back(str6);
          RID.insert(make_pair(str5,v2));
        }
      }
    }
    in.close();
    remove(files[i].c_str());
    vector<string>Uibd;
    for(auto iR = UR.begin(); iR !=UR.end(); iR++)
      Uibd.push_back(iR->first);
    for(int i=0; i<UID.size(); i++){
      if(find(Uibd.begin(), Uibd.end(), UID[i])==Uibd.end()){
        cout<<"individual "<<UID[i]<<" does not exit in the IBDLD file. The program will quit." << endl;
        exit(1);
      }
    }
    for(auto iR = UR.begin(); iR !=UR.end(); iR++){
      int counter =0;
      bool *flagc = new bool[pos9.size()];
      for(int i=0; i<pos9.size(); i++){
        flagc[i] = false;
      }
      vector<vector <double> > v1;
      vector<string> v4;
      v4 = RID[iR->first];
      v1 = iR->second;
      int length;
      length = v1[0].size()-9;
      if(mm==(l-1))
        length = v1[0].size();
#pragma omp parallel for default (none) shared(pos9,i,data2,gpro1,config, \
                                               FU,UR, length,v4,v1,cout,posmap,interval1,freq1, \
                                               genbol1,posibd0,interval,flagc,iR) private(Rdelta, URS, Rbool, Upr) 
      for(int j=0; j<length; j+=9){
        int snp = (j+len1)/9;
        vector<int>v;
        v = interval1[snp];
        for(int k=0; k<v.size(); k++){
          int snp9 = v[k];
          URS.clear();
          Rbool.clear();
          Rdelta.clear();
          vector<double>freq;
          freq= freq1[pos9[snp9]];
          double q=freq[0];
          double p=freq[1];
          int G1, G2, Su1, Su2, S_12;
          int key;
          vector<double>prsu;
          G1= 1;
          G2= 0;
          Su1=8;
          Su2=9;
          S_12=8;
          prsu.push_back(0.0);
          prsu.push_back(.5*q*q);
          prsu.push_back(q*q);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1= 1;
          G2= 0;
          Su1=8;
          Su2=6;
          S_12=5;
          prsu.push_back(0.0);
          prsu.push_back(0.5*q);
          prsu.push_back(q);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1= 1;
          G2= 0;
          Su1=3;
          Su2=4;
          S_12=8;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(q*q);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1= 1;
          G2= 0;
          Su1=3;
          Su2=2;
          S_12=5;
          prsu.push_back(1.0);
          prsu.push_back(1.0);
          prsu.push_back(q);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=1 ;
          G2=2 ;
          Su1=8;
          Su2=9;
          S_12=8;
          prsu.push_back(p*p);
          prsu.push_back(.5*p*p);
          prsu.push_back(0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=1 ;
          G2=2 ;
          Su1=8;
          Su2=6;
          S_12=5;
          prsu.push_back(p);
          prsu.push_back(.5*p);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=1 ;
          G2=2 ;
          Su1=3;
          Su2=4;
          S_12=8;
          prsu.push_back(p*p);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=1 ;
          G2=2 ;
          Su1=3;
          Su2=2;
          S_12=5;
          prsu.push_back(p);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=8;
          Su2=8;
          S_12=9;
          prsu.push_back(q*q);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=5;
          Su2=8;
          S_12=4;
          prsu.push_back(q);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=8;
          Su2=5;
          S_12=6;
          prsu.push_back(q);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=5;
          Su2=5;
          S_12=2;
          prsu.push_back(1.0);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=2 ;
          Su1=8;
          Su2=8;
          S_12=9;
          prsu.push_back(0.0);
          prsu.push_back(0.5*p*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=2 ;
          Su1=5;
          Su2=8;
          S_12=4;
          prsu.push_back(0.0);
          prsu.push_back(0.5*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=2 ;
          Su1=8;
          Su2=5;
          S_12=6;
          prsu.push_back(0.0);
          prsu.push_back(0.5*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=2 ;
          Su1=5;
          Su2=5;
          S_12=2;
          prsu.push_back(0.0);
          prsu.push_back(0.5);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=8;
          Su2=8;
          S_12=9;
          prsu.push_back(0.0);
          prsu.push_back(.5*q*p);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=5;
          Su2=8;
          S_12=4;
          prsu.push_back(0.0);
          prsu.push_back(0.5*p);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=8;
          Su2=5;
          S_12=6;
          prsu.push_back(0.0);
          prsu.push_back(0.5*p);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=5;
          Su2=5;
          S_12=2;
          prsu.push_back(0.0);
          prsu.push_back(0.5);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();



          G1=2 ;
          G2=2 ;
          Su1=8;
          Su2=8;
          S_12=9;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(p*p);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();

          G1=2 ;
          G2=2 ;
          Su1=5;
          Su2=8;
          S_12=4;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(p);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=2 ;
          Su1=8;
          Su2=5;
          S_12=6;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(p);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=2 ;
          Su1=5;
          Su2=5;
          S_12=2;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(1.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();

          G1=0 ;
          G2=0 ;
          Su1=8;
          Su2=8;
          S_12=8;
          prsu.push_back(0.5*q*(q+1));
          prsu.push_back(0.25*q*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=5;
          Su2=8;
          S_12=3;
          prsu.push_back(q);
          prsu.push_back(0.5*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=8;
          Su2=5;
          S_12=5;
          prsu.push_back(q);
          prsu.push_back(0.5*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=2 ;
          Su1=8;
          Su2=8;
          S_12=8;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();

          G1=0 ;
          G2=2 ;
          Su1=5;
          Su2=8;
          S_12=3;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=2 ;
          Su1=8;
          Su2=5;
          S_12=5;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=8;
          Su2=8;
          S_12=8;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=5;
          Su2=8;
          S_12=3;
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=8;
          Su2=5;
          S_12=5;
          prsu.push_back(0);
          prsu.push_back(0);
          prsu.push_back(0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=2 ;
          Su1=8;
          Su2=8;
          S_12=8;
          prsu.push_back(0.0);
          prsu.push_back(0.25*p*p);
          prsu.push_back(0.5*p*(p+1));
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=2 ;
          Su1=5;
          Su2=8;
          S_12=3;
          prsu.push_back(0.0);
          prsu.push_back(0.5*p);
          prsu.push_back(p);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=2 ;
          Su1=8;
          Su2=5;
          S_12=5;
          prsu.push_back(0.0);
          prsu.push_back(0.5*p);
          prsu.push_back(p);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=0 ;
          G2=0 ;
          Su1=8;
          Su2=8;
          S_12=7;
          prsu.push_back(q);
          prsu.push_back(0.5*q);
          prsu.push_back(0.0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();

          G1=0 ;
          G2=2 ;
          Su1=8;
          Su2=8;
          S_12=7;
          prsu.push_back(0);
          prsu.push_back(0);
          prsu.push_back(0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          G1=2 ;
          G2=0 ;
          Su1=8;
          Su2=8;
          S_12=7;
          prsu.push_back(0);
          prsu.push_back(0);
          prsu.push_back(0);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();

          G1=2 ;
          G2=2 ;
          Su1=8;
          Su2=8;
          S_12=7;
          prsu.push_back(0.0);
          prsu.push_back(0.5*p);
          prsu.push_back(p);
          key = G1*1e4+G2*1e3+Su1*1e2+Su2*10+S_12;
          config[key]=prsu;
          prsu.clear();
          double prG1G2S12[9][3][3]={
            {
              {1,    0,    0},
              {0,    1,    0},
              {0,    0,    1}
            },
            {
              {1-p,  0,    p},
              {0,    0,    0},
              {1-p,  0,    p}
            },
            {
              {1,  0,    0},
              {.5, 0,    .5},
              {0,  0,     1}
            },
            {
              {1-p,  0,  p},
              {1-p,  0,  p},
              {1-p,  0,  p}
            },
            {
              {1-p, p,   0},
              {0,   0,   0},
              {0,   1-p, p}
            },
            {
              {(1-p)*(1-p),  p*(1-p),    p*p},
              {(1-p)*(1-p),  p*(1-p),    p*p},
              {(1-p)*(1-p),  2*p*(1-p),  p*p}
            },
            {
              {1, 0,  0},
              {0, 1,  0},
              {0, 0,  1}
            },
            {
              {1-p,      p,    0},
              {.5*(1-p), .5,   .5*p},
              {0,        1-p,  p}
            },
            {
              {(1-p)*(1-p),   2*p*(1-p),  p*p},
              {(1-p)*(1-p),   2*p*(1-p),  p*p},
              {(1-p)*(1-p),   2*p*(1-p),  p*p}
            }
          };
          double prG2S12[9][3]=  {{1-p,         0,         p},
                                  {1-p,         0,         p},
                                  {(1-p)*(1-p),2*p*(1-p), p*p},
                                  {(1-p)*(1-p),2*p*(1-p), p*p},
                                  {1-p,        0,         p},
                                  {1-p,        0,         p},
                                  {(1-p)*(1-p),2*p*(1-p), p*p},
                                  {(1-p)*(1-p),2*p*(1-p), p*p},
                                  {(1-p)*(1-p),2*p*(1-p), p*p}};
          double prsu0;
          double prsu1;
          double prsu2;
          double count=0.0;
          double countcon=0.0;
          double del1, del2, del3, del4, del5, del6, del7, del8, del9;
          for(int ii=0; ii<v1.size();ii++){
            counter++;
            vector<double> v2;
            if(posibd0[0] >pos9[snp9]){
              del1 = v1[ii][0];
              del2 = v1[ii][1];
              del3 = v1[ii][2];
              del4 = v1[ii][3];
              del5 = v1[ii][4];
              del6 = v1[ii][5];
              del7 = v1[ii][6];
              del8 = v1[ii][7];
              del9 = v1[ii][8];
            }else if(posibd0[posibd0.size()-1] <=pos9[snp9]){
              del1 = v1[ii][v1[ii].size()-9];
              del2 = v1[ii][v1[ii].size()-8];
              del3 = v1[ii][v1[ii].size()-7];
              del4 = v1[ii][v1[ii].size()-6];
              del5 = v1[ii][v1[ii].size()-5];
              del6 = v1[ii][v1[ii].size()-4];
              del7 = v1[ii][v1[ii].size()-3];
              del8 = v1[ii][v1[ii].size()-2];
              del9 = v1[ii][v1[ii].size()-1];
            }else{
              del1 = v1[ii][j]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+9]-v1[ii][j])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

              del2 = v1[ii][j+1]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+10]-v1[ii][j+1])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del3 = v1[ii][j+2]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+11]-v1[ii][j+2])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del4 = v1[ii][j+3]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+12]-v1[ii][j+3])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del5 = v1[ii][j+4]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+13]-v1[ii][j+4])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del6 = v1[ii][j+5]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+14]-v1[ii][j+5])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del7 = v1[ii][j+6]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+15]-v1[ii][j+6])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del8 = v1[ii][j+7]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+16]-v1[ii][j+7])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);

              del9 = v1[ii][j+8]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                  (v1[ii][j+17]-v1[ii][j+8])/(posibd0[(interval[snp9]/9)+1]-
                                              posibd0[interval[snp9]/9]);
            }
            v2.push_back(del1);
            v2.push_back(del2);
            v2.push_back(del3);
            v2.push_back(del4);
            v2.push_back(del5);
            v2.push_back(del6);
            v2.push_back(del7);
            v2.push_back(del8);
            v2.push_back(del9);
            Rdelta.insert(make_pair(v4[ii],v2));      
            vector<double> max;
            max_sec(v2,max,9);
            if((max[0]==1 && max[2]==9) || (max[0]==1 && max[2]==4) ||(max[0]==1 && max[2]==6)){}
            else{
              URS.insert(make_pair(v4[ii],max));
              Rbool.insert(make_pair(v4[ii],false));
            }
          }
          for(auto is1 = URS.begin(); is1 !=URS.end(); is1++){
            if(Rbool[is1->first]==false){
              for(auto is2 = is1; is2 !=URS.end(); is2++){
                if(is1->first ==is2->first){}
                else{
                  if(Rbool[is2->first] ==false){
                    double del1, del2, del3, del4, del5, del6, del7, del8, del9;
                    if(posibd0[0] >pos9[snp9]){
                      Key key11={is1->first,is2->first,snp};
                      vector<double>v11;
                      v11 = data2[key11];
                      del1 = v11[0];
                      del2 = v11[1];
                      del3 = v11[2];
                      del4 = v11[3];
                      del5 = v11[4];
                      del6 = v11[5];
                      del7 = v11[6];
                      del8 = v11[7];
                      del9 = v11[8];
                    }else if(posibd0[posibd0.size()-1] <=pos9[snp9]){
                      Key key22={is1->first,is2->first,snp};
                      vector<double>v22;
                      v22 = data2[key22];
                      del1 = v22[0];
                      del2 = v22[1];
                      del3 = v22[2];
                      del4 = v22[3];
                      del5 = v22[4];
                      del6 = v22[5];
                      del7 = v22[6];
                      del8 = v22[7];
                      del9 = v22[8];
                    }else{
                      Key key11={is1->first,is2->first,snp};
                      Key key22={is1->first,is2->first,snp+1};
                      vector<double>v11;
                      vector<double>v22;
                      v11 = data2[key11];
                      v22 = data2[key22];
                      del1 = v11[0]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[0]-v11[0])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);
     
                      del2 = v11[1]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[1]-v11[1])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del3 = v11[2]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[2]-v11[2])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del4 = v11[3]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[3]-v11[3])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del5 = v11[4]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[4]-v11[4])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del6 = v11[5]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[5]-v11[5])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del7 = v11[6]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[6]-v11[6])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del8 = v11[7]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[7]-v11[7])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);

                      del9 = v11[8]+(pos9[snp9]-posibd0[interval[snp9]/9])*
                          (v22[8]-v11[8])/(posibd0[(interval[snp9]/9)+1]-
                                           posibd0[interval[snp9]/9]);
                    }
                    vector<double>v2;
                    v2.push_back(del1);
                    v2.push_back(del2);
                    v2.push_back(del3);
                    v2.push_back(del4);
                    v2.push_back(del5);
                    v2.push_back(del6);
                    v2.push_back(del7);
                    v2.push_back(del8);
                    v2.push_back(del9);
                    vector<double>max;
                    max_sec(v2,max,9);
                    vector<double> v31;
                    v31 = gpro1[make_pair(is1->first,pos9[snp9])];
                    auto l1 = max_element(begin(v31), std::end(v31));
                    int d1 = distance(begin(v31),l1);
                    vector<double> v32;
                    v32 = gpro1[make_pair(is2->first,pos9[snp9])];
                    auto l2 = max_element(begin(v32), std::end(v32));
                    int d2 = distance(begin(v32),l2);
                    int key1 = d1*1e4+d2*1e3+is1->second[2]*1e2+is2->second[2]*10+max[2];
                    int key2 = d2*1e4+d1*1e3+is2->second[2]*1e2+is1->second[2]*10+max[2];
                    auto icon1 = config.find(key1);
                    auto icon2 = config.find(key2);
                    if(icon1 !=config.end()){
                      vector<double>v1;
                      v1 = config[key1];
                      if(key1 == 10898 && flagc[snp9] ==false){
                        flagc[snp9]=true;
                      }else if(key1==10898 && flagc[snp9]==true){
                        v1[0]=v1[0];
                        v1[1]=2*v1[1];
                        v1[2]=v1[2];
                      }
                      if(key1 == 12898 && flagc[snp9] ==false){
                        flagc[snp9]=true;
                      }else if(key1==12898 && flagc[snp9]==true){
                        v1[0]=v1[0];
                        v1[1]=2*v1[1];
                        v1[2]=v1[2];
                      }

                      if(key1 ==20889 || key1 == 2889){
                        v1[0] = 50*v1[0];
                        v1[1] = 50*v1[1];
                        v1[2] = 50*v1[2];
                      }
                      double sum=0.0;
                      for(int i=0; i<2; i++){
                        for(int k=0; k<2; k++){
                          for(int l=0; l<2; l++){
                            sum = sum+is1->second[i]*is2->second[k]*
                                max[l];
                          }
                        }
                      }
                      double vv = is1->second[0]*is2->second[0]*max[0];
                      v1[0] = v1[0]*vv/sum;
                      v1[1] = v1[1]*vv/sum;
                      v1[2] = v1[2]*vv/sum;
                      double sum1=0.0;
                      for(int i=0; i<2; i++){
                        for(int k=0; k<2; k++){
                          for(int l=0; l<2; l++){
                            if(i==0 && k==0 && l==0)
                              continue;
                            sum1= sum1 + prG1G2S12[int(max[2+l])-1][d2][d1]*prG2S12[int(max[2+l])-1][d2]*
                                is1->second[i]*is2->second[k]*max[l]/sum;

                          }
                        }
                      }
                      v1[0] = v1[0]+sum1;
                      v1[1] = v1[1]+sum1;
                      v1[2] = v1[2]+sum1;
                      if(v1[0] ==0)
                        v1[0] = v1[0]+1e-8;
                      if(v1[1] ==0)
                        v1[1] = v1[1]+1e-8;
                      if(v1[2] ==0)
                        v1[2] = v1[2]+1e-8;
                      auto iU = Upr.find(make_pair(iR->first,pos9[snp9]));
                      vector<double>v4;
                      if(iU !=Upr.end()){
                        v4 = iU->second;
                        v4[0] = v4[0]+log(v1[0]);
                        v4[1] = v4[1]+log(v1[1]);
                        v4[2] = v4[2]+log(v1[2]);
                        iU->second = v4;
                      }else{
                        v4.push_back(log(v1[0]));
                        v4.push_back(log(v1[1]));
                        v4.push_back(log(v1[2]));
                        Upr[make_pair(iR->first,pos9[snp9])]= v4;
                      }
                      Rbool[is1->first] = true;
                      Rbool[is2->first] = true;
                      break;
                    }
                    else if(icon2 !=config.end()){
                      vector<double>v1;
                      v1 = config[key2];
                      if(key2 == 10898 && flagc[snp9] ==false){
                        flagc[snp9]=true;
                      }else if(key2==10898 && flagc[snp9]==true){
                        v1[0]=v1[0];
                        v1[1]=2*v1[1];
                        v1[2]=v1[2];
                      }

                      if(key2 == 12898 && flagc[snp9] ==false){
                        flagc[snp9]=true;
                      }else if(key2==12898 && flagc[snp9]==true){
                        v1[0]=v1[0];
                        v1[1]=2*v1[1];
                        v1[2]=v1[2];
                      }

                      if(key2 ==20889 || key2 == 2889){
                        v1[0] = 50*v1[0];
                        v1[1] = 50*v1[1];
                        v1[2] = 50*v1[2];
                      }

                      double sum=0.0;
                      for(int i=0; i<2; i++){
                        for(int k=0; k<2; k++){
                          for(int l=0; l<2; l++){
                            sum = sum+is1->second[i]*is2->second[k]*
                                max[l];

                          }
                        }
                      }
                      double vv = is1->second[0]*is2->second[0]*max[0];
                      v1[0] = v1[0]*vv/sum;
                      v1[1] = v1[1]*vv/sum;
                      v1[2] = v1[2]*vv/sum;
                      double sum1=0.0;
                      for(int i=0; i<2; i++){
                        for(int k=0; k<2; k++){
                          for(int l=0; l<2; l++){
                            if(i==0 && k==0 && l==0)
                              continue;
                            sum1= sum1 + prG1G2S12[int(max[2+l])-1][d1][d2]*prG2S12[int(max[2+l])-1][d1]*
                                is1->second[i]*is2->second[k]*max[l]/sum;
                          }
                        }
                      }
                      v1[0] = v1[0]+sum1;
                      v1[1] = v1[1]+sum1;
                      v1[2] = v1[2]+sum1;
                      if(v1[0] ==0)
                        v1[0] = v1[0]+1e-8;
                      if(v1[1] ==0)
                        v1[1] = v1[1]+1e-8;
                      if(v1[2] ==0)
                        v1[2] = v1[2]+1e-8;
                      auto iU = Upr.find(make_pair(iR->first,pos9[snp9]));
                      vector<double>v4;
                      if(iU !=Upr.end()){
                        v4 = iU->second;
                        v4[0] = v4[0]+log(v1[0]);
                        v4[1] = v4[1]+log(v1[1]);
                        v4[2] = v4[2]+log(v1[2]);
                        iU->second = v4;
                      }else{
                        v4.push_back(log(v1[0]));
                        v4.push_back(log(v1[1]));
                        v4.push_back(log(v1[2]));
                        Upr[make_pair(iR->first,pos9[snp9])]= v4;
                      }
                      Rbool[is1->first] = true;
                      Rbool[is2->first] = true;
                      break;

                 
                    }
                  }
                }
              }
            }

            if(Rbool[is1->first]==false){
              auto l1 = max_element(begin(gpro1[make_pair(is1->first,pos9[snp9])]), std::end(gpro1[make_pair(is1->first,pos9[snp9])]));
              int d1 = distance(begin(gpro1[make_pair(is1->first,pos9[snp9])]),l1);

              vector<double> v;
              v = Rdelta[is1->first];
              double f= v[0]+v[1]+v[4]+v[5];
              auto iFU = FU.find(make_pair(iR->first,snp9));
              if(iFU !=FU.end()){
                iFU->second = iFU->second+f;
              }else{
                FU[make_pair(iR->first,snp9)]=f;
              }  
              Upr2(Rdelta[is1->first], Upr, p, d1, f, pos9[snp9], iR->first);
              Rbool[is1->first] = true;
            }
          }
        }
      }
      delete [] flagc;
      if(counter < SID.size())
        cout<<" Only "<<counter<<" individuals used to impute the individual "<<iR->first<<endl; 
    }
    len1 = len1 + len2;
  }

  ofstream output1(myfile8);
  double counterr=0.0;
  double counter22=0.0;
  if(flag){
    ifstream ifile22("ppf");
    while(getline(ifile22, line1)){
      stringstream s3(line1);
      s3 >> str1 >> str2 >> val1 >> str3 >> str4;
      vector<double> v1;
      double val22;
      while(s3 >> val22){
        v1.push_back(val22);
      }
      output1<<"---"<<" "<<chro<<":"<<val1<<":"<<allele[val1][0]<< 
          ":"<<allele[val1][1]<<" "<<val1<<" "<<allele[val1][0]<<" "<<
          allele[val1][1]<<" ";
      for(int i =0; i<v1.size(); i=i+3){
        double pru02 = v1[i];
        double  pru12 = v1[i+1];
        double  pru22 = v1[i+2];
        vector<double>v2;
        v2 = Upr[make_pair(UID[i/3],val1)];
        if(v2.size() ==0){
          cout<< "Upr1 is zero "<<v2.size()<<endl;
          exit(1);
        }
        vector<double>vv;
        vv.push_back(pru02);
        vv.push_back(pru12);
        vv.push_back(pru22);
        auto l4 = max_element(begin(vv), std::end(vv));
        double d4 = distance(begin(vv),l4);
        if(pru02*exp(v2[0])==0 && pru12*exp(v2[1])==0 && pru22*exp(v2[2])==0){
          if(vv[0]==1){
            pru02 = vv[0]-0.01;
          }else if(vv[0]==0.0){
            pru02 = vv[0]+0.005;
          }else{pru02 = vv[0];}
          if(vv[1] ==1){
            pru12 = vv[1]-0.01;
          }else if(vv[1] ==0.0){
            pru12 = vv[1]+0.005;
          }else{pru12 =vv[1];}
          if(vv[2] ==1){
            pru22 = vv[2]-0.01;
          }else if(vv[2] ==0.0){
            pru22 = vv[2]+0.005;
          }else{pru22 =vv[2];}
        }
        vector<double>freq;
        freq = freq1[val1];
        double p = freq[1];
        double fu =  Uf[make_pair(UID[i/3],val1)];
        double pru0= fu*(1-p)+(1-fu)*pow((1-p),2.0);
        double pru1= (1-fu)*2*p*(1-p);
        double pru2= fu*p+(1-fu)*pow(p,2.0);
        double sum = 0.0;
        sum = exp(v2[0])*pru0 + exp(v2[1])*pru1 + exp(v2[2])*pru2;
        double pr0 = exp(v2[0])*pru0/sum;
        double pr1 = exp(v2[1])*pru1/sum;
        double pr2 = exp(v2[2])*pru2/sum;
        vector<double> pr;
        pr.push_back(pr0);
        pr.push_back(pr1);
        pr.push_back(pr2);
        double sum2 = 0.0;
        sum2 = exp(v2[0])*pru02 + exp(v2[1])*pru12 + exp(v2[2])*pru22;
        double pr02 = exp(v2[0])*pru02/sum2;
        double pr12 = exp(v2[1])*pru12/sum2;
        double pr22 = exp(v2[2])*pru22/sum2;
        vector<double> prr;
        prr.push_back(pr02);
        prr.push_back(pr12);
        prr.push_back(pr22);
        auto l2 = max_element(begin(prr), std::end(prr));
        double d2 = distance(begin(prr),l2);
        auto l = max_element(begin(pr), std::end(pr));
        double d = distance(begin(pr),l);
        output1<<pr02<<" "<<pr12<<" "<<pr22<<" ";
        if(*l2 >=.9)
          counter22++;
      }
      output1<<endl;
    }
    cout<<"call con rate prior i2 "<<counter22/Upr.size()<<endl;
  }else{
    for(int& pos:pos9){
      output1<<"---"<<" "<<chro<<":"<<pos<<":"<<allele[pos][0]<<
          ":"<<allele[pos][1]<<" "<<pos<<" "<<allele[pos][0]<<" "<<
          allele[pos][1]<<" ";
      for(string& str:UID){
        vector<double>v1;
        v1 = Upr[make_pair(str,pos)];
        if(v1.size() ==0){
          cout<< "Upr is zero "<<v1.size()<<endl;
          exit(1);
        }
        vector<double>freq;
        freq = freq1[pos];
        double p = freq[1];
        double fu =  Uf[make_pair(str,pos)];
        double pru0= fu*(1-p)+(1-fu)*pow((1-p),2.0);
        double pru1= (1-fu)*2*p*(1-p);
        double pru2= fu*p+(1-fu)*pow(p,2.0);
        double sum = 0.0;
        sum = exp(v1[0])*pru0 + exp(v1[1])*pru1 + exp(v1[2])*pru2;
        double pr0 = exp(v1[0])*pru0/sum;
        double pr1 = exp(v1[1])*pru1/sum;
        double pr2 = exp(v1[2])*pru2/sum;
        vector<double> pr;
        pr.push_back(pr0);
        pr.push_back(pr1);
        pr.push_back(pr2);
        auto l = max_element(begin(pr), std::end(pr));
        double d = distance(begin(pr),l);
        output1<<pr0<<" "<<pr1<<" "<<pr2<<" ";
        if(*l >=.9)   
          counterr++;
      }
      output1 <<endl;
    }
  }

  //     allinea_stop_sampling();
}

void get_delta(vector<double>& delta1, vector<double>& delta2, vector<double>& delta3, vector<double>& delta4, vector<double>& delta5, vector<double>& delta6,
               vector<double>& delta7, vector<double>& delta8, vector<double>& delta9,
               vector<double>& v, int Ns){
  int l1=0;
  for(int i=0; i<Ns; i++){
    delta1.push_back(v[l1]);
    delta2.push_back(v[1+l1]);
    delta3.push_back(v[2+l1]);
    delta4.push_back(v[3+l1]);
    delta5.push_back(v[4+l1]);
    delta6.push_back(v[5+l1]);
    delta7.push_back(v[6+l1]);
    delta8.push_back(v[7+l1]);
    delta9.push_back(v[8+l1]);
    l1 +=9;
  }
}
void get_W(vector<double>&W0, vector<double>& W1, vector<double>& W2,double delta1,
           double delta2, double delta3, double delta4, double delta5, double delta6,
           double delta7, double delta8, double delta9, double p){
  W0.push_back((1-p)*delta1);
  W0.push_back((1-p)*delta2);
  W0.push_back(pow((1-p),2.0)*delta3);
  W0.push_back(pow((1-p),2.0)*delta4);
  W0.push_back((1-p)*delta5);
  W0.push_back((1-p)*delta6);
  W0.push_back(pow((1-p),2.0)*delta7);
  W0.push_back(pow((1-p),2.0)*delta8);
  W0.push_back(pow((1-p),2.0)*delta9);
  W1.push_back(0.0*delta1);
  W1.push_back(0.0*delta2);
  W1.push_back(2*p*(1-p)*delta3);
  W1.push_back(2*p*(1-p)*delta4);
  W1.push_back(0.0*delta5);
  W1.push_back(0.0*delta6);
  W1.push_back(2*p*(1-p)*delta7);
  W1.push_back(2*p*(1-p)*delta8);
  W1.push_back(2*p*(1-p)*delta9);
  W2.push_back(p*delta1);
  W2.push_back(p*delta2);
  W2.push_back(pow(p,2.0)*delta3);
  W2.push_back(pow(p,2.0)*delta4);
  W2.push_back(p*delta5);
  W2.push_back(p*delta6);
  W2.push_back(pow(p,2.0)*delta7);
  W2.push_back(pow(p,2.0)*delta8);
  W2.push_back(pow(p,2.0)*delta9);
}

void vec_mat_mult(double e[],double A[][9],int m, int n,double eA[]){
  for(int i=0; i<n; i++){
    eA[i]=0.0;
    for(int j=0 ; j<m; j++){
      eA[i]+= e[j]*A[j][i];
    }
  }
}

void Upr2(vector<double>& v, unordered_map<US_pair, vector<double>, hash_name1>& Upr,double p, double gen, double fu, int snp, string U){
  vector<double>W0;
  vector<double>W1;
  vector<double>W2;
  double epsilon = 0.005;
  get_W(W0, W1, W2, v[0], v[1], v[2], v[3], v[4], v[5], v[6],v[7], v[8], p);

  double A0[3][9] = {{1, 1-p, 1, 1-p, 1-p, pow((1-p),2.0), 1, 1-p, pow((1-p),2.0)},
                     {0, 0  , 0, 0  , p  , p*(1-p)       , 0, p  , 2*p*(1-p)},
                     {0, p  , 0, p  , 0  , pow(p,2)      , 0, 0  , pow(p,2.0)}};
  double A1[3][9] ={{0, 0, .5, 1-p, 0, pow((1-p),2.0), 0, .5*(1-p), pow((1-p),2.0)},
                    {1, 0, 0 ,  0 , 0, p*(1-p)       ,1 , .5      , 2*p*(1-p)     },
                    {0, 0, .5,  p , 0, pow(p,2.0),    0, .5*p,     pow(p,2.0)     }};
  double A2[3][9] ={{0, 1-p, 0, 1-p, 0 , pow((1-p),2.0), 0, 0  ,  pow((1-p),2.0)},
                    {0, 0  , 0, 0 , 1-p, 2*p*(1-p),      0, 1-p, 2*p*(1-p)      },
                    {1, p  , 1, p , p, pow(p,2.0) ,      1,  p , pow(p,2.0)     }};
  double e0[3]= {pow((1-epsilon),2.0), epsilon*(1-epsilon), pow(epsilon,2)};
  double e1[3]= {2*epsilon*(1-epsilon), 1-2*epsilon*(1-epsilon), 2*epsilon*(1-epsilon)};
  double e2[3]= {pow(epsilon,2.0), epsilon*(1-epsilon), pow((1-epsilon),2.0)};
  int m=3;
  int n=9;
  double pru0= fu*(1-p)+(1-fu)*pow((1-p),2.0);
  double pru1= (1-fu)*2*p*(1-p);
  double pru2= fu*p+(1-fu)*pow(p,2.0);
  if(pru0 == 0){
    pru0 = pru0+1e-2;
  }
  if(pru1 == 0){
    pru1 = pru1+1e-2;
  }
  if(pru2 == 0){
    pru2 = pru2+1e-2;
  }
  if(gen==0){
    double e0A0[9];
    double e0A1[9];
    double e0A2[9];
    vec_mat_mult(e0,A0,m,n,e0A0);
    vec_mat_mult(e0,A1,m,n,e0A1);
    vec_mat_mult(e0,A2,m,n,e0A2);
    double eAW0=0.0, eAW1=0.0, eAW2=0.0;
    for(int i=0; i<9; i++){
      eAW0 += e0A0[i]*W0[i];
      eAW1 += e0A1[i]*W1[i];
      eAW2 += e0A2[i]*W2[i];
    }
    double prsu0 = log((1/(pru0))*eAW0);
    double prsu1 = log((1/(pru1))*eAW1);
    double prsu2 = log((1/(pru2))*eAW2);
    auto iU = Upr.find(make_pair(U,snp));
    vector<double>v6;
    if(iU !=Upr.end()){
      v6 = iU->second;
      v6[0] = v6[0]+prsu0;
      v6[1] = v6[1]+prsu1;
      v6[2] = v6[2]+prsu2;
      iU->second = v6;
    }else{
      v6.push_back(prsu0);
      v6.push_back(prsu1);
      v6.push_back(prsu2);
      Upr[make_pair(U,snp)]= v6;
    }
  }else if(gen==1){
    double e1A0[9];
    double e1A1[9];
    double e1A2[9];
    vec_mat_mult(e1,A0,m,n,e1A0);
    vec_mat_mult(e1,A1,m,n,e1A1);
    vec_mat_mult(e1,A2,m,n,e1A2);
    double eAW0=0.0, eAW1=0.0, eAW2=0.0;
    for(int i=0; i<9; i++){
      eAW0 += e1A0[i]*W0[i];
      eAW1 += e1A1[i]*W1[i];
      eAW2 += e1A2[i]*W2[i];
    }
    double prsu0 = log((1/(pru0))*eAW0);
    double prsu1 = log((1/(pru1))*eAW1);
    double prsu2 = log((1/(pru2))*eAW2);
    auto iU = Upr.find(make_pair(U,snp));
    vector<double>v6;
    if(iU !=Upr.end()){
      v6 = iU->second;
      v6[0] = v6[0]+prsu0;
      v6[1] = v6[1]+prsu1;
      v6[2] = v6[2]+prsu2;
      iU->second = v6;
    }else{
      v6.push_back(prsu0);
      v6.push_back(prsu1);
      v6.push_back(prsu2);
      Upr[make_pair(U,snp)]= v6;
    }
  }else if(gen==2){
    double e2A0[9];
    double e2A1[9];
    double e2A2[9];
    vec_mat_mult(e2,A0,m,n,e2A0);
    vec_mat_mult(e2,A1,m,n,e2A1);
    vec_mat_mult(e2,A2,m,n,e2A2);
    double eAW0=0.0, eAW1=0.0, eAW2=0.0;
    for(int i=0; i<9; i++){
      eAW0 += e2A0[i]*W0[i];
      eAW1 += e2A1[i]*W1[i];
      eAW2 += e2A2[i]*W2[i];
    }
    double prsu0 = log((1/(pru0))*eAW0);
    double prsu1 = log((1/(pru1))*eAW1);
    double prsu2 = log((1/(pru2))*eAW2);
    auto iU = Upr.find(make_pair(U,snp));
    vector<double>v6;
    if(iU !=Upr.end()){
      v6 = iU->second;
      v6[0] = v6[0]+prsu0;
      v6[1] = v6[1]+prsu1;
      v6[2] = v6[2]+prsu2;
      iU->second = v6;
    }else{
      v6.push_back(prsu0);
      v6.push_back(prsu1);
      v6.push_back(prsu2);
      Upr[make_pair(U,snp)]= v6;
    }
  }
}

void max_sec(vector<double>& v3, vector<double>& max, int N){
  double maxi, secmax;
  int maxindex, secmaxindex;
  if(v3[0]>v3[1]){
    maxi = v3[0];
    maxindex =0;
    secmax= v3[1];
    secmaxindex = 1;
  }else{
    maxi = v3[1];
    maxindex = 1;
    secmax = v3[0];
    secmaxindex = 0;
  }
  for(int i=2; i<N; i++){
    if(v3[i] >maxi){
      secmax = maxi;
      maxi = v3[i];
      secmaxindex = maxindex;
      maxindex = i;
    }else if(v3[i]>secmax && v3[i] !=maxi){
      secmax =v3[i];
      secmaxindex = i;
    }
  }
  maxindex++;
  secmaxindex++;
  max.push_back(maxi);
  max.push_back(secmax);
  max.push_back(maxindex);
  max.push_back(secmaxindex);
}
void prGUS(int S, double gen, unordered_map<US_pair, vector<double>, hash_name1>& Upr, int snp, string U, double prG1G2S12[][3][3]){
  if(gen ==1000){}
  else{
    double prsu0 = prG1G2S12[S-1][0][int(gen)];
    double prsu1 = prG1G2S12[S-1][1][int(gen)];
    double prsu2 = prG1G2S12[S-1][2][int(gen)];
    if(prsu0 ==0)
      prsu0 = prsu0 +1e-8;
    if(prsu1 ==0)
      prsu1 = prsu1 +1e-8;
    if(prsu2 ==0)
      prsu2 = prsu2 +1e-8;

    vector<double>v4;
    auto iU = Upr.find(make_pair(U,snp));
    if(iU !=Upr.end()){
      v4 = iU->second;
      v4[0] = v4[0]+log(prsu0);
      v4[1] = v4[1]+log(prsu1);
      v4[2] = v4[2]+log(prsu2);
      iU->second = v4;
    }else{
      v4.push_back(log(prsu0));
      v4.push_back(log(prsu1));
      v4.push_back(log(prsu2));
      Upr[make_pair(U,snp)]= v4;
    }
  }

}




//    ifile110.clear();
//    ifile110.seekg(0, ios::beg);
