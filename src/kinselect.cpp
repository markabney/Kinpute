/*
 * Kinselect is a program to select individuals to sequence that optimizes
 * later genotype imputation.
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
#include <limits>
//#include <time.h>
//#include <string>
#include <cstring>

using namespace std;

typedef pair<string,string> ID_pair;

struct hash_name {
  size_t operator()(const ID_pair &name ) const
  {
    return hash<string>()(name.first) ^ hash<string>()(name.second);
  }
};

double Lmin(vector<int>& S, vector<string>& pop,
            unordered_map<string,int>& ID, double** w, int N)
{
  double Lnew = 0.0;
  double prod;
  for(string& p : pop){
    prod = 1;
    for(int k=0; k < S.size(); k++){
      if(ID[p] == S[k]){
        prod = 0.0;
        break;
      }else{
        prod = prod * (1 - w[S[k]][ID[p]]);
      }
    }
    Lnew = Lnew + prod;
  }
  return Lnew;
}

int main (int argc,char* argv[]) {
  string myfile1;
  string myfile2;
  string myfile3;
  string samplefile;
  int M;
  bool verbose = false;
  
  if(argc < 6){
    cout << "Usage is: " << argv[0] << " -i <infile> -or <outfile1> -ou <outfile2> "
         << " -n N -s <samplefile> [-v]" << endl;
    exit(0);
  }else{
    for(int i=1; i < argc; i++){
      //if(i+1 != argc){
      if(strcmp(argv[i], "-i")==0){
        myfile1 = argv[i+1];
      }else if(strcmp(argv[i], "-or")==0){
        myfile2 = argv[i+1];
      }else if(strcmp(argv[i], "-ou")==0){
        myfile3 = argv[i+1];
      } else if (strcmp(argv[i], "-n")==0){
        M =stoi(argv[i+1]);
      } else if (strcmp(argv[i], "-s")==0){
        samplefile = argv[i+1];
      } else if (strcmp(argv[i], "-v")==0) {
        verbose = true;
      }
      //}
    }
  }
  
  unordered_map <string, int> ID ;
  unordered_map <string, int> :: iterator found;
  
  int N = 0;
  vector<string> pop;
  ifstream sampfl(samplefile);
  if (!sampfl) {
    cerr << "Could not open " << samplefile << "\n";
    exit(1);
  }
  string inLine;
  while (getline(sampfl, inLine)) {
    stringstream sline(inLine);
    string name;
    sline >> name;

    found = ID.find(name);
    if (found != ID.end()) {
      if (verbose) {
        cerr << "Warning: Individual " << name << " appears multiple times in file "
             << samplefile << ". Ignoring the repeated name.\n";
      }        
      continue;
    } 
    pop.push_back(name);
    ID.insert(make_pair(name, N));
    ++N;
  }
  
  double Lold, Lnew;
  int vertex;
  double **w;
  w = new double *[N];
  for (int i=0; i<N; i++)
    w[i] = new double[N];

  // unordered_map <string, int> :: iterator it ;
  unordered_map<ID_pair,double,hash_name> Kinship;

  ifstream ifile1 (myfile1);
  string name1, name2;
  double kincoef;
  int k=0;
  while(ifile1 >> name1 >> name2 >> kincoef) {
    if (ID.find(name1) != ID.end() && ID.find(name2) != ID.end()) {
      ID_pair p1(name1, name2);
      Kinship.insert(make_pair(p1, kincoef));
      ID_pair p2(name2, name1);
      Kinship.insert(make_pair(p2, kincoef));      
    }

  }

  for(string& p1:pop){
    for(string& p2:pop){
      if (p1 == p2) continue;
      if (Kinship.find(ID_pair(p1, p2)) != Kinship.end()) {
        w[ID[p1]][ID[p2]] = Kinship[ID_pair(p1, p2)];  
      } else {
        if (verbose) {
          cerr << "Warning: No kinship coefficient found for pair " << p1
               << " " << p2 << " . Assuming a value of 0.\n";  
        }        
        w[ID[p1]][ID[p2]] = 0;
      }
    }
  }

  vector <int> S;
  Lold = std::numeric_limits<double>::max();

  for(int l=0; l<M; l++){
    for(string& p:pop){
      if(std::find(S.begin(), S.end(), ID[p]) != S.end())
      {}
      else {
        S.push_back(ID[p]);
        Lnew = Lmin(S, pop, ID, w, N);
        S.pop_back();
        if(Lnew < Lold){
          Lold = Lnew;
          vertex = ID[p];
        }
      }
    }
    S.push_back(vertex);
  }

  ofstream outputr;
  outputr.open(myfile2);

  ofstream outputu;
  outputu.open(myfile3);
  for(int i=0; i<S.size(); i++){
    for(auto it=ID.begin(); it !=ID.end(); it++){
      if(it->second ==S[i])
        outputr <<it->first<<endl;
    }
  }

  vector<int> U;
  for(int i=0; i<N; i++){
    if(std::find(S.begin(), S.end(), i) !=S.end()){}
    else{ U.push_back(i);}
  }

  for(int i=0; i<U.size(); i++){
    for(auto it=ID.begin(); it !=ID.end(); it++){
      if(it->second ==U[i])
        outputu <<it->first<<endl;
    }
  }





}

