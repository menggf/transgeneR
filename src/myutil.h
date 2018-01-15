/*
  myutil.cpp @ transgeneR
  Provide the shared function and format the output
  
  Copyright (c) 2017- Guofeng Meng
  BT science biotechnology
  
  EMAIL: menggf@gmail.com
*/
#include<iostream>
#include<ostream>
#include<string>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<cmath>
#include<sstream>
#include<map>
//#include <Rcpp.h>
#include <stdlib.h>
#include<string.h>


using namespace std;

vector <string> split(string s, const string delim="\t") {
    vector<string> result;
    size_t pos=0;
    while(1){
        size_t temp=s.find(delim, pos);
        if(temp==string::npos){
            result.push_back(s.substr(pos,1000));
            return result;
        }
        if(temp==pos)
            result.push_back("");
        else
            result.push_back(s.substr(pos,temp-pos));
        pos=temp+1;
    }
    return result;
}
string joinstr(vector <string> v, const string delim="\t")
{
	stringstream ss;
	for(size_t i = 0; i < v.size(); ++i)
	{
  	if(i != 0)
   		ss << delim;
  		ss << v[i];
	}

	string s = ss.str();
	return s;
}
bool dfd(map<string, string> a, string b){
	map <string,string>::iterator it = a.find(b);
  	if (it != a.end())
  		return true;
  	else
  		return false;
}

bool dfd(map<int, string> a, int b){
	map <int,string>::iterator it = a.find(b);
  	if (it != a.end())
  		return true;
  	else
  		return false;
}

bool dfd(map<int, int> a, int b){
	map <int,int>::iterator it = a.find(b);
  	if (it != a.end())
  		return true;
  	else
  		return false;
}
bool dfd(map<char, int> a, char b){
	map <char,int>::iterator it = a.find(b);
  	if (it != a.end())
  		return true;
  	else
  		return false;
}



int toi(string x){
	return atoi(x.c_str());
}
template <typename T>
string to_string ( T Number )
{
     std::ostringstream ss;
     ss << Number;
     return ss.str();
}
  
int gts(int n){
	if(n>=4096) 
		return 4096;
	if(n>=2048)
		return 2048;
	if(n>=1024)
		return 1024;
	if(n>=512)
		return 512;
	if(n>=256)
		return 256;
	if(n>=128) 
		return 128;
	if(n>=64)	
		return 64;
	if(n>=32)
		return 32;
	if(n>=16)
		return 16;
	if(n>=8)
		return 8;
	if(n>=4) 
		return 4;
	if(n>=2) 
		return 2;
	if(n>=1)
		return 1;
	return 0;
}

map <int, int> flags(int n){
	map <int,int> res;
	res[1]=0;
	res[2]=0;
	res[4]=0;
	res[8]=0;
	res[16]=0;
	res[32]=0;
	res[64]=0;
	res[128]=0;
	res[256]=0;
	res[512]=0;
	res[1024]=0;
	res[2048]=0;
	res[4096]=0;
	while(n > 0){
		int x=gts(n);
		res[x]=1;
		n=n-x;
	}
	return res;
}
template <typename T>
T mymin(T a, T b){
	if(a > b)
		return b;
	else
		return a;
}
template <typename T>
T mymax(T a, T b){
	if(a > b)
		return a;
	else
		return b;
}
int issplit(int from1, int to1,int from2, int to2, int total=150){
	if(from1 <= from2)
		 if(to1 < from2 + 10 && to1-from1 + to2-from2 + 2 > 0.7 * float(total))
			return 1;
	if(from1 > from2)
		if(to2 < from1 + 10 && to1-from1 + to2-from2 + 2 > 0.7 * float(total))
			return -1;
	return 0;
}
int isoverlap(int from1, int to1, int from2, int to2){
	if(from1 <= from2)
		 if(to1 > from2 && to1-from2 > mymin(10.0, 0.8* mymin(to1-from1, to2-from2)))
			return 1;
	if(from1 > from2)
		 if(to2 > from1 && to2-from1 > mymin(10.0, 0.8* mymin(to1-from1, to2- from2)))
			return 1;
	return 0;
}


//judeg the match parts of reads are in the homo region
int ishomo(string chr0, int from0, int to0, vector <string> chr, vector <int> from1, vector <int> to1, vector <int> from2, vector <int> to2){ 
	int tag=0;
	int num=from1.size();
	for(int i=0; i < num; i++){
		if(chr[i] != chr0)
			continue;
		if(from0 > from1[i]-30 && from0 < to1[i] && to0 > from1[i] && to0 < to1[i]+ 30){
			return 1;
		}
	}
	return tag;
}

map <char, int> regx(string x){
	int num=x.length();
	char * y = new char [x.length()+1];
  	strcpy (y, x.c_str());
	map <char, int> res;
	int l=0;
	int z[5];
	z[0]=0;
	z[1]=0;
	z[2]=0;
	z[3]=0;
	z[4]=0;
	int wh=0;
	for(int i=0;i<num;i++){
		if(y[i] <=57 && y[i] >=48){
			z[l]=y[i]-48;
			l++;
		}
		else{
			int sum=0;
			for(int j=0;j < l;j++){
				sum=sum+ z[j] *(pow(10, (l-j-1)));
				z[j]=0;
			}
			l=0;
			wh++;
			if(y[i]=='S' && wh==1){
				res['1']=sum;
			}
			else if(y[i]=='S')
				res['2']=sum;
			else{
				if(dfd(res, y[i]))
					res[y[i]]+=sum;
				else
					res[y[i]]=sum;
			}
		}
	}
	delete y;
	//cout<<x<<endl;
	//for (map<char,int>::iterator it=res.begin(); it!=res.end(); ++it)
    //	cout << it->first << " => " << it->second << '\n';
	return res;
}

int len_cigar(string cigar){
	map <char, int> reg=regx(cigar);
	int sum=0;
	map<char,int>::iterator it;
	for (it=reg.begin(); it!=reg.end(); ++it)
    	if(it->first == 'M')
    		sum+= it->second;
    	if(it->first == 'I')
    		sum+= it->second;
    	if(it->first == 'N')
    		sum+=it->second;
    	if(it->first == '=')
    		sum+=it->second;
    	if(it->first == 'X')
    		sum+=it->second;
	return sum;
}

void extract_soft(string cigar, int *ll){
	map <char, int> reg=regx(cigar);
	if(!dfd(reg,'1'))
		ll[0]=0;
	else
		ll[0]=reg['1'];
	if(!dfd(reg,'2'))
		ll[2]=0;
	else
		ll[2]=reg['2'];
	ll[1]=len_cigar(cigar);
}
int scores(char *idstr){
	int sc=-1000;
	char * pch = strstr(idstr,"AS:i:");
	if(pch ==NULL)
		return sc;
	int w=0;
	char tmp[10];
	for(w=0;w<10;w++){
		if(pch[w+5] == ' ')
			break;
		tmp[w]=pch[w+5];
		
	}
	sc=atoi(tmp);
	return sc;
}
int scores2(char *idstr){
	int sc=-1000;
	char * pch = strstr(idstr,"XS:i:");
	if(pch ==NULL)
		return sc;
	int w=0;
	char tmp[10];
	for(w=0;w<10;w++){
		if(pch[w+5] == ' ')
			break;
		tmp[w]=pch[w+5];
		
	}
	sc=atoi(tmp);
	return sc;
}
int xm(char *idstr){
	int sc=100;
	char * pch = strstr(idstr,"XM:i:");
	if(pch ==NULL)
		return sc;
	int w=0;
	char tmp[10];
	for(w=0;w<10;w++){
		if(pch[w+5] == ' ')
			break;
		tmp[w]=pch[w+5];
		
	}
	sc=atoi(tmp);
	return sc;
}


string revff(string ff){
	if(ff=="ff")
		return "rr";
	if(ff=="fr")
		return "fr";
	if(ff=="rf")
		return "rf";
	if(ff=="rr")
		return "ff";
	return ff;
}

bool msc( string chr1, string from,  string chr2, string to, string ff){ // change the order or not?
	if(chr1==chr2){
		if(toi(from)> toi(to))
			return(true);
		else
			return(false);
	}
	else if(chr1=="insert" && chr2!="insert" )
		return(true);
	else if(chr1!="insert" && chr2=="insert" )
		return(false);
	else if(chr1 !="insert" && chr2!="insert"){
		if(chr1 < chr2)
			return(false);
		else
			return(true);
	}
	return(false);
}

string ms(bool is, string name, string chr1, string from, string gap, string chr2, string to, string ff,string sc, string assign){
	if(is)
		return("insert\t"+name+"\t"+chr2+"\t"+to+"\t"+gap+"\t"+chr1+"\t"+from+"\t"+revff(ff) +"\t"+sc + "\t" + assign);
	else
		return("insert\t"+name+"\t"+chr1+"\t"+from+"\t"+gap+"\t"+chr2+"\t"+to+"\t"+ff+"\t"+sc + "\t" + assign);
}

string mc(bool is, string name, string ff, string chr1, string from1, string to1,string chr2, string from2, string to2,string chr3, string from3, string to3,string chr4, string from4, string to4){
	if(is)
		return("split\t"+name +"\t"+revff(ff)+"\t"+chr4+"\t"+from4+"\t"+to4+"\t"+chr3+"\t"+from3+"\t"+to3+"\t"+chr2+"\t"+from2+"\t"+to2+"\t"+chr1+"\t"+from1+"\t"+to1);
	else
		return("split\t"+name +"\t"+ff+"\t"+chr1+"\t"+from1+"\t"+to1+"\t"+chr2+"\t"+from2+"\t"+to2+"\t"+chr3+"\t"+from3+"\t"+to3+"\t"+chr4+"\t"+from4+"\t"+to4);
}

vector <string> readreads(char * idstr){
	vector <string> res=split(idstr,"\t");
	map <int,int> tg=flags(toi(res[1]));
	string x="0";
	if(tg[64]==1)
		x="1";
	if(tg[128]==1)
		x="2";
	vector <string> rr;
  int sc= scores(idstr);
  int sc2= scores2(idstr);
	if(tg[4]==1){
		rr.push_back(x);
		rr.push_back(string(""));
		return rr;
	}
	int ll[3];
    extract_soft(res[5],ll);
        
    string strand="1";
	if(tg[16]==1)
		strand= "-1";
    string cord="N";
    if(res[8] != "0")
    	cord="Y";
    if(cord == "Y")
		if(abs(toi(res[3])-toi(res[7])) > 1000)
			cord="N" ;
	
	int nm = xm(idstr);
	
	string temp=res[2]+"\t"+res[3]+"\t"+strand+"\t"+cord+"\t"+res[5]+"\t"+res[8]+"\t"+res[9]+"\t"+to_string(sc)+"\t"+to_string(ll[0])+"\t"+to_string(ll[1])+"\t"+to_string(ll[2])+"\t"+to_string(nm) +"\t"+to_string(sc-sc2);
	rr.push_back(x);
	rr.push_back(temp);
	return(rr);
}





