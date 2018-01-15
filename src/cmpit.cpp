/*
  reassignreads.cpp @ transgeneR
  Reads assembly using the second-round alignment results.

  Copyright (c) 2017- Guofeng Meng
  BT science biotechnology

  EMAIL: menggf@gmail.com
*/
#include<iostream>
#include<string>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<cmath>
#include<sstream>
#include<map>
#include <Rcpp.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

using namespace std;

int findfirst(int i, int c1, int c2,  Rcpp::IntegerVector arr){
  int from=0;
  int to=c2-1;
  while(from < to){
    int p=int((from + to)/2);
    if(arr(3+i) > int(arr(3+ c1 + c2 +p)))
      from=p;
    else if(arr(3+i) < int(arr(3+ c1 +p)))
      to=p;
    else
      break;
    if(arr(3+i) <= int(arr(3+ c1 + c2 +p)) && arr(3+i) >= int(arr(3+ c1 +p)))
      return from;
    if(to - from <= 3000 )
      return from;
  }
  return 0;
}
// [[Rcpp::export]]
Rcpp::IntegerVector cmpitcpp(Rcpp::IntegerVector arr){
  int const cores=int(arr(0));
	int const c1=int(arr(1));
	int const c2=int(arr(2));
	Rcpp::IntegerVector cc(c1, 0);
//  int nProcessors=cores;
//	omp_set_num_threads(nProcessors);
//#pragma omp parallel for
	for(int i=0; i < c1; i++ ){
     int j=findfirst(i, c1, c2,arr);
		for(; j< c2;j++){
			if(int(arr(3+i)) > int(arr(3+c1+c2+j)))
				continue;
			if(int(arr(3+i)) <= int(arr(3+c1+c2+j))){
				if(int(arr(3+i)) >= int(arr(3+c1+j))){
//#pragma omp critical
//{
					cc[i]++;
//}
        }
			}
			if(int(arr(3+i)) < int(arr(3+c1+j)) - 100)
				break;
		}
	}
	return(cc);
}

int maxfreq(vector <int> x){
  sort(x.begin(), x.end());
	int tt=-1, cc=1;
	int wh_max=-1,  cc_max=-1;
	int n=x.size();
	while(tt < n-1){
     tt++;

		if(x[tt] == -1)
			continue;
		if( tt != n -1 && x[tt] == x[tt+1] ){
			cc++;
		}
		else{
			if(cc > cc_max){
				wh_max=x[tt];
				cc_max=cc;
			}
			cc=0;
		}
	}
	return wh_max;
}

// [[Rcpp::export]]
Rcpp::List clustercpp(Rcpp::IntegerVector arr){
  const int n=arr.size() -1 ;
	Rcpp::IntegerVector cl(n,  0);
	vector <int> ct;
	int const rg=int(arr(0));
	vector <int> x;
	for(int i=1;i<= n;i++ )
		x.push_back(arr(i));


	int zz=0;

	while(1){
		int wh=maxfreq(x);
     if(zz > 100)
       break;
     if(wh==-1)
       break;
		zz++;
		ct.push_back(wh);
		for(int i=0; i< n; i++){
			if(x[i] >= wh -rg && x[i] <= wh + rg && x[i] != -1 ){
				cl[i]=zz;
				x[i]=-1;
			}
		}
	}
	Rcpp::IntegerVector center(ct.size());
	for(int i=0;i<ct.size(); i++)
		center[i]=ct[i];
	return Rcpp::List::create(Rcpp::Named("cl") = cl,Rcpp::Named("center") = center);
}
