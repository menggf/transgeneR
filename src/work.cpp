/*
  assignreads.cpp @ transgeneR
  Reads assignment to the first-round alignment results.

  Copyright (c) 2017- Guofeng Meng
  BT science biotechnology

  EMAIL: menggf@gmail.com
*/



// [[Rcpp::depends(BH)]]
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
#include "myutil.h"
#include <string.h>
#include <omp.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


using namespace std;
vector <string> dohomo(string dir, int win, int step){
  vector <string> all;
  char idstr[1000];
  map <string, string> seq;
  fstream fseq((string(dir)+"/temp_homo.fa").c_str(),fstream::in);
  while(fseq.getline(idstr,1000)){
    int idstr_len=sizeof(idstr)/sizeof(char);
    char id[1000];
    memcpy(id, &idstr[1], idstr_len-2 );
    fseq.getline(idstr,1000);
    seq[string(id)]=idstr;
  }
  fseq.close();

  vector <string> chr;
  vector <int> std, f1, f2, t1, t2;
  fstream fhomo((string(dir)+"/temp_homo.sam").c_str(),fstream::in);
  while(fhomo.getline(idstr,1000)){
    if(idstr[0] =='@')
      continue;
    vector <string> res=split(idstr,"\t");
    map <int,int> tg=flags(toi(res[1]));
    if(tg[4]==1)
      continue;

    int ll[3];
    extract_soft(res[5],ll);
    int strand=0;
    if(seq[res[0]] == res[9])
      strand= 1;
    else
      strand= -1;
    chr.push_back(res[2]);
    std.push_back(strand);
    if(strand==1){
      f1.push_back(toi(res[3]));
      t1.push_back(toi(res[3])+ll[1]);
      f2.push_back(toi(res[0])+ll[0]+1);
      t2.push_back(toi(res[0])+ ll[0]+ll[1]);
    }
    else{
      t1.push_back(toi(res[3]));
      f1.push_back(toi(res[3])+ll[1]);
      f2.push_back(toi(res[0])+ ll[2]+1);
      t2.push_back(toi(res[0])+ ll[2]+ll[1]);
    }
  }
  fhomo.close();

  int num=std.size();
  int istrand=0;
  string ichr="";
  int ifrom1=0;
  int ifrom2=0;
  int ito1=0;
  int ito2=0;
  for(int i=0;i<num;i++){
    if(ichr==""){
      istrand=std[i];
      ichr=chr[i];
      ifrom1=f1[i];
      ifrom2=f2[i];
      ito1=t1[i];
      ito2=t2[i];
    }
    else{
      if(istrand == std[i] && ichr == chr[i] && abs(abs(ito2-f2[i]) - win + step) < 10 && abs(abs(ito1-f1[i])-win+step)< 10){
        ito1=t1[i];
        ito2=t2[i];
      }
      else{
        if(istrand==1)
          all.push_back(ichr+"\t"+to_string(istrand)+"\t"+to_string(ifrom1)+"\t"+to_string(ito1)+"\t"+to_string(ifrom2)+"\t"+to_string(ito2));
        else
          all.push_back(ichr+"\t"+to_string(istrand)+"\t"+to_string(ito1)+"\t"+to_string(ifrom1)+"\t"+to_string(ifrom2)+"\t"+to_string(ito2));

        istrand=std[i];
        ichr=chr[i];
        ifrom1=f1[i];
        ifrom2=f2[i];
        ito1=t1[i];
        ito2=t2[i];
      }
    }
  }
  return all;
}


// [[Rcpp::export]]
void findhomo(Rcpp::String dir, Rcpp::IntegerVector para){
  vector <string> output=dohomo(dir, int(para(0)),int(para(1)));
  fstream fout1((string(dir)+"/homo.txt").c_str(),fstream::out);
  int num=output.size();
  for(int i=0;i<num;i++)
    fout1<<output[i]<<endl;
  fout1.close();
}


void work(vector <string> &all, vector <string> &used_ids, vector <string> &read1, vector <string> &read2, vector <string> &read3,vector <string> &read4, vector <string> &homo_chr, vector <int> &homo_from1, vector <int> &homo_to1, vector <int> &homo_from2, vector <int> &homo_to2, int min_frac, int strict=1, int tail=10){
	int num_ids=used_ids.size();
    for(int xx=0; xx < num_ids; xx++){
    	int id=xx;
    	//if(used_ids[id]!="E00477:175:H7V7WCCXx:6:2220:12540:33129")
    	//	continue;
    	string assign="unknown";
    	if(read1[id]=="" && read2[id]=="" && read3[id]=="" && read4[id]=="")
    		continue;
    	if(read1[id]!="" && read2[id]=="" && read3[id] !="" && read4[id]=="")
    		continue;
    	if(read1[id]=="" && read2[id] !="" && read3[id] =="" && read4[id] !="")
    		continue;
    	if(read1[id]!="" && read2[id] !="" && read3[id] !="" && read4[id] !=""){ //1234
    		vector <string> ss1=split(read1[id],"\t");
    		vector <string> ss2=split(read2[id],"\t");
    		vector <string> ss3=split(read3[id],"\t");
    		vector <string> ss4=split(read4[id],"\t");
    		int dup1, dup2, splt1,splt2;
    		if(ss1[6] == ss3[6]){ // the same direction
    			// the same query sequence?
				dup1= isoverlap(toi(ss1[8]), toi(ss1[8])+toi(ss1[9]), toi(ss3[8]), toi(ss3[8])+toi(ss3[9]));
				// mapped split reads
				splt1= issplit(toi(ss1[8]), toi(ss1[8])+toi(ss1[9]), toi(ss3[8]), toi(ss3[8])+toi(ss3[9]), ss1[6].length());
    		}
    		else{ // the different direction
				dup1=   isoverlap(toi(ss1[8]), toi(ss1[8])+toi(ss1[9]), toi(ss3[10]), toi(ss3[10])+toi(ss3[9]));
				splt1=  issplit(toi(ss1[8]), toi(ss1[8])+toi(ss1[9]), toi(ss3[10]), toi(ss3[10])+toi(ss3[9]), ss1[6].length());
    		}
    		if(ss2[6] == ss4[6]){ // the same direction
				dup2=   isoverlap(toi(ss2[8]), toi(ss2[9])+toi(ss2[8]),toi(ss4[8]), toi(ss4[8]) + toi(ss4[9]));
				splt2=  issplit(toi(ss2[8]), toi(ss2[9])+toi(ss2[8]), toi(ss4[8]), toi(ss4[8])+toi(ss4[9]), ss2[6].length());
    		}
    		else{ // the different direction
				dup2=   isoverlap(toi(ss2[8]), toi(ss2[9])+toi(ss2[8]),toi(ss4[10]), toi(ss4[10])+toi(ss4[9]));
				splt2=  issplit(toi(ss2[8]), toi(ss2[9])+toi(ss2[8]),toi(ss4[10]), toi(ss4[10])+toi(ss4[9]), ss2[6].length());
    		}

			int splt3= issplit(toi(ss1[8]), toi(ss1[8])+toi(ss1[9]), toi(ss2[8]), toi(ss2[8])+toi(ss2[9]), ss1[6].length());
			int splt4= issplit(toi(ss3[8]), toi(ss3[8])+toi(ss3[9]), toi(ss4[8]), toi(ss4[8])+toi(ss4[9]), ss1[6].length());
			int dup3=   isoverlap(toi(ss1[8]), toi(ss1[8])+toi(ss1[9]), toi(ss2[8]), toi(ss2[8])+toi(ss2[10]));
			int dup4=   isoverlap(toi(ss3[8]), toi(ss3[8])+toi(ss3[9]), toi(ss4[8]), toi(ss4[8])+toi(ss4[10]));
			//if(dup3  && splt3 )
			//	continue;

			int ishomo1=ishomo(ss1[0], toi(ss1[1]), toi(ss1[1])+ toi(ss1[9]), homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
			int ishomo2=ishomo(ss2[0], toi(ss2[1]), toi(ss2[1])+ toi(ss2[9]), homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
			//cout<<used_ids[id]<<"\t"<<ishomo1<<" "<<ishomo2<<" "<< dup1 <<" "<< dup2<<" "<< splt1<<" "<< splt2<<endl;
			if((ishomo1 || ishomo2) && (toi(ss1[8]) > tail || toi(ss2[8]) > tail || toi(ss1[10]) > tail || toi(ss2[10]) > tail)){
				read1[id]="";
				read2[id]="";
			}
			else if(ss1[3] == "Y" && ss3[3] == "Y"){
				int sc1=toi(ss1[7]) + toi(ss2[7]);
				int sc2=toi(ss3[7]) + toi(ss4[7]);
				if(dup1  && dup2){
					if( sc1 > sc2 + 10 ){
						read3[id]="";
						read4[id]="";
					}
					else if(sc2 > sc1 + 10){
						read1[id]="";
						read2[id]="";
					}
					else if(toi(ss1[8]) > tail || toi(ss1[10]) > tail || toi(ss2[8]) > tail || toi(ss2[10]) > tail){
						read1[id]="";
						read2[id]="";
					}
					else{
						string xx2="assignboth\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[1]+"\t"+to_string(toi(ss1[1]) + toi(ss1[9]));
						string xx3="assignboth\t"+used_ids[id]+"\t"+ss2[0]+"\t"+ss2[1]+"\t"+to_string(toi(ss2[1]) + toi(ss2[9]));
						string xx4="assignboth\t"+used_ids[id]+"\t"+ss3[0]+"\t"+ss3[1]+"\t"+to_string(toi(ss3[1]) + toi(ss3[9]));
						string xx5="assignboth\t"+used_ids[id]+"\t"+ss4[0]+"\t"+ss4[1]+"\t"+to_string(toi(ss4[1]) + toi(ss4[9]));
						all.push_back(xx2);
						all.push_back(xx3);
						all.push_back(xx4);
						all.push_back(xx5);
						continue;
					}
				}
				else if(splt1 && splt2){
					string ff="ff";
					if(toi(ss1[5]) < 0){
						ss1=split(read2[id],"\t");
						ss2=split(read1[id],"\t");
						ss3=split(read4[id],"\t");
						ss4=split(read3[id],"\t");
					}
					if(splt1==-1 && splt2==-1){
						if(ss1[6] == ss3[6])
							ff="ff";
						else if(ss1[6] != ss3[6])
							ff="rf";
						else{
						  Rcpp::Rcout<<used_ids[id]<<" "<<ss1[2]<<" "<<ss3[2]<<endl;
						  Rcpp::Rcout<<"4Error?"<<endl;
						}

						int gap=toi(ss1[8]) - toi(ss3[9]) - toi(ss3[8]);
						int to1=toi(ss3[1]) + toi(ss3[9]);
						if(ff == "fr" || ff == "rf"){
							gap=toi(ss1[8]) - toi(ss3[9])- toi(ss3[10]);
							to1=toi(ss3[1]);
						}
						bool is =msc(ss3[0],to_string(to1),ss1[0],ss1[1],ff);
						string xx1=ms(is, used_ids[id],ss3[0],to_string(to1),to_string(gap),ss1[0],ss1[1],ff,to_string(sqrt(toi(ss1[7]) * toi(ss3[7]))), "splitreads");
						all.push_back(xx1);
						string xx3;
						if(toi(ss1[1]) + toi(ss1[9]) >= toi(ss2[1]) && ss1[3] == "Y")
							xx3=mc(is,used_ids[id],ff,"","","",ss3[0],ss3[1],to_string(toi(ss3[1]) + toi(ss3[9])),ss1[0],ss1[1],to_string(toi(ss2[1]) + toi(ss2[9])),"","","");
						else
							xx3=mc(is,used_ids[id],ff,"","","",ss3[0],ss3[1],to_string(toi(ss3[1]) + toi(ss3[9])),ss1[0],ss1[1],to_string(toi(ss1[1]) + toi(ss1[9])),ss2[0],ss2[1],to_string(toi(ss2[1]) + toi(ss2[9])));
						all.push_back(xx3);
						continue;
					}
					else if(splt1 == 1 && splt2 == 1){
						if(ss1[6] == ss3[6])
							ff="ff";
						else if(ss1[6] != ss3[6])
							ff="fr";
						else{
						  Rcpp::Rcout<<used_ids[id]<<" "<<ss1[2]<<" "<<ss3[2]<<endl;
						  Rcpp::Rcout<<"1Error?"<<endl;
						}
						int gap=toi(ss2[10]) - toi(ss4[9]) - toi(ss4[10]);
						int from2=toi(ss2[1]) + toi(ss2[9]);
						int from4=toi(ss4[1]);
						if(ff == "fr" || ff == "rf"){
							gap=toi(ss2[10])- toi(ss4[9]) - toi(ss4[8]);
							from4=toi(ss4[1]) + toi(ss4[9]);
						}
						bool is=msc(ss2[0],to_string(from2),ss4[0],to_string(from4),ff);
						string xx1=ms(is, used_ids[id],ss2[0],to_string(from2),to_string(gap),ss4[0],to_string(from4),ff,to_string(sqrt(toi(ss2[7])*toi(ss4[7]))), "splitreads");
						all.push_back(xx1);
						string xx2;
						if(toi(ss3[1]) + toi(ss3[9]) >= toi(ss4[1]) && ss2[3] == "Y")
							xx2=mc(is, used_ids[id],ff,"", "", "", ss1[0],ss1[1],to_string(toi(ss1[1]) + toi(ss1[9])),ss3[0],ss3[1],to_string(toi(ss4[1]) + toi(ss4[9])),"","","");
						else
							xx2=mc(is, used_ids[id],ff,"", "", "", ss1[0],ss1[1],to_string(toi(ss1[1]) + toi(ss1[9])),ss3[0],ss3[1],to_string(toi(ss3[1]) + toi(ss3[9])),ss4[0],ss4[1],to_string(toi(ss4[1]) + toi(ss4[9])));
						all.push_back(xx2);
						continue;
					}
					else{
						all.push_back("warning\t"+used_ids[id]+": split direction looks strange:1234!");
						continue;
					}
					all.push_back("assign\t"+used_ids[id]+"\tsplitreads");
					all.push_back("assign2\t"+used_ids[id]+"\tsplitreads");
					assign="splitreads";
					continue;
				}
				else{
					all.push_back("warning\t"+used_ids[id]+": unknown how to assign!");
					continue;
				}
			}
			else if(ss1[3] == "Y"  && ss3[3] == "N"){
				if(splt1 && !splt2)
					read4[id]="";
				else if(!splt1 && splt2)
					read3[id]="";
				else if(!splt1 && !splt2){
					read3[id]="";
					read4[id]="";
				}
				else if(splt1 && splt2){
					if(rand() ==1){
						read3[id]="";
					}
					else{
						read4[id]="";
					}
				}
				else{
					all.push_back("warning\t"+used_ids[id]+": looks strange: NY!");
					continue;
				}
			}
			else if(ss1[3] == "N"  && ss3[3] == "Y"){
				if(ishomo1 || ishomo2){
					if(ishomo1){
						read1[id]="";
					}
					if(ishomo2){
						read2[id]="";
					}
				}
				else if(splt1 && !splt2){
					read1[id]="";
				}
				else if(!splt1 && splt2){
					read1[id]="";
				}
				else if(ss1[7] < ss3[7] && ss2[7] < ss4[7]){
					read1[id]="";
					read1[id]="";
				}
				else{
					all.push_back("warning\t"+used_ids[id]+": looks strange: NY!");
					continue;
				}
			}
			else if(ss1[3] == "N"  && ss3[3]  == "N"){  // both not cordinated mapped
				if(toi(ss1[7]) > toi(ss3[7])){
					read3[id]="";
				}
				else{
					read1[id]="";
				}
				if(toi(ss2[7]) > toi(ss4[7])){
					read4[id]="";
				}
				else{
					read2[id]="";
				}
			}
			else{
				all.push_back("warning\t"+used_ids[id]+": did not find the assignment: 1234!");
			}
    	}
		if(read1[id]!="" && read2[id]!=""){ // change the direction for first time
			vector <string> ss1=split(read1[id],"\t");
 		  vector <string> ss2=split(read2[id],"\t");
			if(ss1[3] == "Y"){
				if(toi(ss1[5]) < 0){
					ss1[5]="1";
					ss2[5]="-1";
					read1[id]=joinstr(ss2);
					read2[id]=joinstr(ss1);
					if(read3[id]!="" && read4[id]==""){
						read4[id]=read3[id];
						read3[id]="";
					}
					else if(read3[id]=="" && read4[id]!=""){
						read3[id]=read4[id];
						read4[id]="";
					}
					else if(read3[id]!="" && read4[id]!=""){
						all.push_back("warning\t"+used_ids[id]+": something is wrong: change direction?");
						continue;
					}
				}
			}
		}
		else if(read3[id]!="" && read4[id]!=""){
			vector <string> ss3=split(read3[id],"\t");
    		vector <string> ss4=split(read4[id],"\t");
			if(ss3[3] == "Y"){
				if(toi(ss3[5]) < 0){
					ss3[5]="1";
					ss4[5]="-1";
					read3[id]=joinstr(ss4);
					read4[id]=joinstr(ss3);
					if(read1[id]!=""){
						read2[id]=read1[id];
						read1[id]="";
					}
					else if(read2[id]!=""){
						read1[id]=read2[id];
						read2[id]="";
					}
				}
			}
		}
		if((read1[id]!="" && read3[id]!="" && read2[id]!="" && read4[id]=="") || (read1[id]!="" && read3[id]!="" && read2[id]=="" && read4[id]!="")){ // not 4 or 2
			vector <string> ss1, ss2,ss3;
			string miss="4";
			if(read4[id]==""){
				ss1=split(read1[id],"\t");
				ss2=split(read2[id],"\t");
				ss3=split(read3[id],"\t");
			}
			else if(read2[id]==""){
				miss="2";
				ss1=split(read3[id],"\t");
				ss2=split(read4[id],"\t");
				ss3=split(read1[id],"\t");
			}
			int dup, splt;
			string ff="ff";

			if(ss1[6] == ss3[6]){ // the same direction
				dup= isoverlap(toi(ss1[8])+1, toi(ss1[8])+toi(ss1[9]), toi(ss3[8])+1, toi(ss3[8])+toi(ss3[9])); // the same query sequence?
				splt=issplit(toi(ss1[8])+1, toi(ss1[8])+toi(ss1[9]), toi(ss3[8])+1, toi(ss3[8])+toi(ss3[9]), ss1[6].length()); // mapped split reads
			}
			else if(ss1[6] != ss3[6] ){ // the different direction
				dup= isoverlap(toi(ss1[8])+1, toi(ss1[8])+toi(ss1[9]), toi(ss3[10])+1, toi(ss3[10])+toi(ss3[9]));
				splt=issplit(toi(ss1[8])+1, toi(ss1[8])+toi(ss1[9]), toi(ss3[10])+1, toi(ss3[10])+toi(ss3[9]), ss1[6].length());
			}
			else{
				all.push_back("warning\t"+used_ids[id]+": Has no mapped strand annotation:N"+miss+"!");
				continue;
			}
			//cout<<used_ids[id]<<" N24 "<<ff<<" "<< dup <<" "<< splt<<endl;
			int ishomo3=ishomo(ss3[0], toi(ss3[1]), toi(ss3[1])+toi(ss3[9]), homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
			if(miss=="2" && ishomo3 && (toi(ss3[8]) > tail || toi(ss3[10]) > tail )){
				read1[id]="";
			}
			else if(dup){
				if(ss1[3] == "Y" || toi(ss1[7]) >= toi(ss3[7]))
					read3[id]="";
				else if(toi(ss1[7]) < toi(ss3[7]))
					read1[id]="";
				else{
					string xx3="assignboth\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[1]+"\t"+to_string(toi(ss1[1]) + toi(ss1[9]));
					string xx4="assignboth\t"+used_ids[id]+"\t"+ss2[0]+"\t"+ss2[1]+"\t"+to_string(toi(ss2[1]) + toi(ss2[9]));
					string xx5="assignboth\t"+used_ids[id]+"\t"+ss3[0]+"\t"+ss3[1]+"\t"+to_string(toi(ss3[1]) + toi(ss3[9]));
					all.push_back(xx3);
					all.push_back(xx4);
					//all.push_back(xx5);
					//all.push_back("assign\t"+used_ids[id]+"\t"+"assignboth");
					continue;
				}
			}
			else if(splt){
				if(ss1[6] != ss3[6]){
					if(splt==-1)
						ff="rf";
					else{
						if(ss1[3] == "Y" ){
							all.push_back("warning\t"+used_ids[id]+": Not used for unknown splits:N" + miss + "!");
							continue;
						}
						else{
							ff="fr";
						}
					}
				}
				if(ss1[3] == "N" && strict){
					all.push_back("warning\t"+used_ids[id]+": look strange in cordination annotation: N"+miss+"!");
					continue;
				}
				if(splt ==1){
					all.push_back("warning\t"+used_ids[id]+": This looks impossible:N"+ miss+"!");
					continue;
				}
				int from=toi(ss3[1]) + toi(ss3[9]);
				int to=toi(ss1[1]);
				int gap=toi(ss1[8]) - toi(ss3[8]) - toi(ss3[9]);
				if(ff == "fr" || ff == "rf"){
					from=toi(ss3[1]);
					gap=toi(ss1[8])-toi(ss3[10])-toi(ss3[9]);
				}
				assign="splitread1";
				if(miss=="2")
					assign="splitread3";
				bool is=msc(ss3[0],to_string(from), ss1[0],to_string(to),ff);
				string xx1=ms(is, used_ids[id],ss3[0],to_string(from),to_string(gap),ss1[0],to_string(to),ff,to_string(sqrt(toi(ss1[7])*toi(ss3[7]))), assign);
				string xx3;
				if(toi(ss1[1]) + toi(ss1[9]) >= toi(ss2[1]) && ss1[3] == "Y")
					xx3=mc(is, used_ids[id],ff,"","","",ss3[0],ss3[1],to_string(toi(ss3[1]) + toi(ss3[9])),ss1[0],ss1[1],to_string(toi(ss2[1]) + toi(ss2[9])),"","","");
				else
					xx3=mc(is, used_ids[id],ff,"","","",ss3[0],ss3[1],to_string(toi(ss3[1]) + toi(ss3[9])),ss1[0],ss1[1],to_string(toi(ss1[1]) + toi(ss1[9])),ss2[0],ss2[1],to_string(toi(ss2[1]) + toi(ss2[9])));

				all.push_back(xx1);
				all.push_back(xx3);

				all.push_back("assign\t"+used_ids[id]+"\t"+assign);
        all.push_back("assign2\t"+used_ids[id]+"\t"+assign);
				continue;
			}
			else{
				all.push_back("warning\t"+used_ids[id]+": Not used for unknown splits:N" + miss + "!");
				continue;
			}
		}
		if( (read1[id]!="" && read3[id]=="" && read2[id]!="" && read4[id]!="") || (read1[id]=="" && read3[id]!="" && read2[id]!="" && read4[id]!="")){ //not 3 or 1
			vector <string> ss1, ss2, ss4;
			string miss="3";
			if(read3[id]==""){
				ss1=split(read1[id],"\t");
				ss2=split(read2[id],"\t");
				ss4=split(read4[id],"\t");
			}
			else if(read1[id]==""){
				miss="1";
				ss1=split(read3[id],"\t");
				ss2=split(read4[id],"\t");
				ss4=split(read2[id],"\t");
			}
			string ff="ff";

			int dup, splt;
			if(ss2[6] == ss4[6]){ // the same direction
				dup= isoverlap(toi(ss2[8])+1, toi(ss2[8])+toi(ss2[9]), toi(ss4[8])+1, toi(ss4[8])+toi(ss4[9])); // the same query sequence?
				splt=  issplit(toi(ss2[8])+1, toi(ss2[8])+toi(ss2[9]), toi(ss4[8])+1, toi(ss4[8])+toi(ss4[9]),ss2[6].length()); // mapped split reads
			}
			else if(ss2[6] != ss4[6]){ // the different direction
				dup= isoverlap(toi(ss2[8])+1, toi(ss2[8])+toi(ss2[9]), toi(ss4[10])+1, toi(ss4[10])+toi(ss4[9])); //the same query sequence?
				splt=  issplit(toi(ss2[8])+1, toi(ss2[8])+toi(ss2[9]), toi(ss4[10])+1, toi(ss4[10])+toi(ss4[9]),ss2[6].length()); //mapped split reads
			}
			//cout<<used_ids[id]<<" N13 "<<ff<<" "<< dup <<" "<< splt<<endl;
			int ishomo2=ishomo(ss4[0], toi(ss4[1]), toi(ss4[1])+toi(ss4[9]), homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
			if(miss=="1" && ishomo2 && (toi(ss4[8]) > tail || toi(ss4[10]) > tail )){
				read2[id]="";
			}
			else if(dup){
				if(ss2[3] == "Y" || toi(ss2[7]) >= toi(ss4[7]))
					read4[id]="";
				else if(toi(ss2[7]) < toi(ss4[7]))
					read2[id]="";
				else{
					string xx3="assignboth\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[1]+"\t"+to_string(toi(ss1[1]) + toi(ss1[9]));
					string xx4="assignboth\t"+used_ids[id]+"\t"+ss2[0]+"\t"+ss2[1]+"\t"+to_string(toi(ss2[1]) + toi(ss2[9]));
					string xx5="assignboth\t"+used_ids[id]+"\t"+ss4[0]+"\t"+ss4[1]+"\t"+to_string(toi(ss4[1]) + toi(ss4[9]));
					all.push_back(xx3);
					all.push_back(xx4);
					//all.push_back(xx5);
					continue;
				}
			}
			else if(splt){
				if(ss2[6] != ss4[6]){
					if(splt == 1)
						ff="fr";
					else
					{
						if(ss2[3] == "Y"){
							all.push_back("warning\t"+used_ids[id]+": Not used for unknown splits:N" + miss + "!");
							continue;
						}
						else
							ff="rf";
					}
				}
				if(ss2[3] == "N" && strict){
					all.push_back("warning\t"+used_ids[id]+ ": look strange in codination annotation:N" + miss + "!");
					continue;
				}
				if(splt == -1){
					all.push_back("warning\t"+used_ids[id]+": This looks impossible:N" + miss + "!");
					continue;
				}
				int from=toi(ss2[1]) + toi(ss2[9]);
				int to=toi(ss4[1]);
				int gap=toi(ss2[10]) - toi(ss4[9])-toi(ss4[10]);
				if(ff == "fr" || ff == "rf"){
					gap=toi(ss2[10]) - toi(ss4[9])- toi(ss4[8]);
					to=toi(ss4[1])+ toi(ss4[9]);
				}
				assign="splitread2";
				if(miss=="1")
					assign="splitread4";
				bool is=msc(ss2[0],to_string(from),ss4[0],to_string(to),ff);
				string xx1=ms(is, used_ids[id],ss2[0],to_string(from),to_string(gap),ss4[0],to_string(to),ff,to_string(sqrt(toi(ss2[7])*toi(ss4[7]))), assign);
				string xx3;
				if(toi(ss1[1]) + toi(ss1[9]) >= toi(ss2[1]))
					xx3=mc(is,used_ids[id],ff,"","","",ss1[0], ss1[1], to_string(toi(ss2[1]) + toi(ss2[9])),ss4[0],ss4[1],to_string(toi(ss4[1]) + toi(ss4[9])), "","","");
				else
					xx3=mc(is, used_ids[id],ff,ss1[0],ss1[1],to_string(toi(ss1[1]) + toi(ss1[9])),ss2[0],ss2[1],to_string(toi(ss2[1]) + toi(ss2[9])),ss4[0],ss4[1],to_string(toi(ss4[1]) + toi(ss4[9])), "","","");

				all.push_back(xx1);
				all.push_back(xx3);
				all.push_back("assign\t"+used_ids[id]+"\t"+assign);
        all.push_back("assign2\t"+used_ids[id]+"\t"+assign);
				continue;
			}
			else{
				all.push_back("warning\t"+used_ids[id]+": Not used for unknown splits:N" + miss + "!");
				continue;
			}
		}
		if( (read1[id]!="" && read3[id]=="" && read2[id]=="" && read4[id]!="") || (read1[id]=="" && read3[id]!="" && read2[id]!="" and read4[id]=="")){ //  1 and 4 or 2 and 3
			vector <string> ss1,ss4;
			int tag=1;
			if(read2[id]==""){
				ss1=split(read1[id],"\t");
				ss4=split(read4[id],"\t");
				assign="splitmid14";
			}
			if(read1[id]==""){
				ss1=split(read3[id],"\t");
				ss4=split(read2[id],"\t");
				assign="splitmid23";
			}
			all.push_back("assign\t"+used_ids[id]+"\t"+assign);
			all.push_back("assign2\t"+used_ids[id]+"\t"+assign);
			if(toi(ss1[8]) > min_frac){
				all.push_back("seq\t>"+used_ids[id]+"_read1left");
				string s=ss1[6].substr(0, toi(ss1[8]));
				all.push_back("seq\t"+s);
				all.push_back("read1left\t"+used_ids[id]+"\t"+s);
				tag=0;
			}
			if(toi(ss1[10]) > min_frac){
				all.push_back("seq\t>"+used_ids[id]+"_read1right");
				string s=ss1[6].substr(toi(ss1[8])+toi(ss1[9]), toi(ss1[10]));
				all.push_back("seq\t"+s);
				all.push_back("read1right\t"+used_ids[id]+"\t"+s);
				tag=0;
			}
			if(toi(ss4[8]) > min_frac){
				all.push_back("seq\t>"+used_ids[id]+"_read2left");
				string s=ss4[6].substr(0, toi(ss4[8]));
				all.push_back("seq\t"+s);
				all.push_back("read2left\t"+used_ids[id]+"\t"+s);
				tag=0;
			}
			if(toi(ss4[10]) > min_frac){
				all.push_back("seq\t>"+used_ids[id]+"_read2right");
				string s=ss4[6].substr(toi(ss4[8])+toi(ss4[9]), toi(ss4[10]));
				all.push_back("seq\t"+s);
				all.push_back("read2right\t"+used_ids[id]+"\t"+s);
				tag=0;
			}
			if(tag){
				all.push_back("range\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[1]+"\t"+to_string(toi(ss1[1])+toi(ss1[9]))+"\t"+ss4[0]+"\t"+ss4[1]+"\t"+to_string(toi(ss4[1])+toi(ss4[9])));
				continue;
			}
			int to1=toi(ss1[1]) + toi(ss1[9]);
			int to4=toi(ss4[1]) + toi(ss4[9]);
			string xx="mark\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[2]+"\t"+ss1[8]+"\t"+ss1[1]+"\t"+to_string(to1)+"\t"+ss1[10]+"\t"+ss1[7]+"\t"+ss4[0]+"\t"+ss4[2]+"\t"+ss4[8]+"\t"+ss4[1]+"\t"+to_string(to4)+"\t"+ss4[10]+"\t"+ss4[7];
			all.push_back(xx);
			continue;
		}
		if((read3[id]=="" && read4[id]=="" && read1[id]!="" && read2[id]!="") || (read1[id]=="" && read2[id]=="" && read3[id]!="" && read4[id]!="")){ // only in genome or insert
			string where="genome";
			vector <string> ss1, ss2;
			if(read1[id]!=""){
				ss1=split(read1[id],"\t");
				ss2=split(read2[id],"\t");
			}
			if(read3[id]!=""){
				where="insert";
				ss1=split(read3[id],"\t");
				ss2=split(read4[id],"\t");
			}
			if(ss1[3] == "Y"){
				assign=where+"concord";
				int tag=1;
				if(toi(ss1[8]) > mymin(min_frac, 40)){
					all.push_back("seq\t>"+used_ids[id]+"_read1left");
					string s=ss1[6].substr(0, toi(ss1[8]));
					all.push_back("seq\t"+s);
					all.push_back("read1left\t"+used_ids[id]+"\t"+s);
					tag=0;
				}
				if(toi(ss2[10]) > mymin(min_frac, 40)){
					all.push_back("seq\t>"+used_ids[id]+"_read2right");
					string s=ss2[6].substr(toi(ss2[8])+toi(ss2[9]), toi(ss2[10]));
					all.push_back("seq\t"+s);
					all.push_back("read2right\t"+used_ids[id]+"\t"+s);
					tag=0;
				}
				if(toi(ss1[8]) > tail & tag == 1){
					if(toi(ss1[1])+toi(ss1[9]) >= toi(ss2[1])){
						all.push_back("nonsplitleft\t"+used_ids[id]+"\tff\t"+ss1[0]+"\t"+ss1[1]+"\t"+ to_string(toi(ss2[1])+toi(ss2[9])));
					}
					else{
						all.push_back("nonsplitleft\t"+used_ids[id]+"\tff\t"+ss1[0]+"\t"+ss1[1]+"\t"+ to_string(toi(ss1[1])+toi(ss1[9])));
                        all.push_back("nonsplitleft\t"+used_ids[id]+"\tff\t"+ss2[0]+"\t"+ss2[1]+"\t"+ to_string(toi(ss2[1])+toi(ss2[9])));
					}
					all.push_back("assign\t"+used_ids[id]+"\t"+assign);
					continue;
				}
				if(toi(ss2[10]) > tail & tag == 1){
					assign=where+"concordmapped";
					if(toi(ss1[1])+toi(ss1[9]) >= toi(ss2[1])){
						all.push_back("nonsplitright\t"+used_ids[id]+"\tff\t"+ss1[0]+"\t"+ss1[1]+"\t"+ to_string(toi(ss2[1])+toi(ss2[9])));
					}
					else{
						all.push_back("nonsplitright\t"+used_ids[id]+"\tff\t"+ss1[0]+"\t"+ss1[1]+"\t"+ to_string(toi(ss1[1])+toi(ss1[9])));
                        all.push_back("nonsplitright\t"+used_ids[id]+"\tff\t"+ss2[0]+"\t"+ss2[1]+"\t"+ to_string(toi(ss2[1])+toi(ss2[9])));
					}
					all.push_back("assign\t"+used_ids[id]+"\t"+assign);
					continue;
				}
				if(tag==1){
					assign=where+"concordmapped";
					if(toi(ss1[1]) > toi(ss2[1])){
						all.push_back("warning\t"+used_ids[id]+ ": reads location seems strange, skipped!");
						continue;
					}
					if(toi(ss1[1])+toi(ss1[9]) >= toi(ss2[1])){
						all.push_back("nonsplit\t"+used_ids[id]+"\tff\t"+ss1[0]+"\t"+ss1[1]+"\t"+ to_string(toi(ss2[1])+toi(ss2[9])));
					}
					else{
						all.push_back("nonsplit\t"+used_ids[id]+"\tff\t"+ss1[0]+"\t"+ss1[1]+"\t"+ to_string(toi(ss1[1])+toi(ss1[9])));
                        all.push_back("nonsplit\t"+used_ids[id]+"\tff\t"+ss2[0]+"\t"+ss2[1]+"\t"+ to_string(toi(ss2[1])+toi(ss2[9])));
					}
					all.push_back("assign\t"+used_ids[id]+"\t"+assign);
					continue;
				}
				all.push_back("assign\t"+used_ids[id]+"\t"+assign);
				all.push_back("assign2\t"+used_ids[id]+"\t"+assign);
				if(toi(ss1[1]) + toi(ss1[9])  > toi(ss2[1])){
					int to=toi(ss2[1]) + toi(ss2[9]);
					string xx="mark\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[2]+"\t"+ss1[8]+"\t"+ss1[1]+"\t"+to_string(to)+"\t"+ss2[10]+"\t"+ss2[7];
					all.push_back(xx);
				}
				else{
					int to2=toi(ss2[1]) + toi(ss2[9]);
					int to1=toi(ss1[1]) + toi(ss1[9]);
					string xx="mark\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[2]+"\t"+ss1[8]+"\t"+ss1[1]+"\t"+to_string(to1)+"\t"+ss1[10]+"\t"+ss1[7]+"\t"+ss2[0]+"\t"+ss2[2]+"\t"+ss2[8]+"\t"+ss2[1]+"\t"+to_string(to2)+"\t"+ss2[10]+"\t"+ss2[7];
					all.push_back(xx);
				}
			}
			else{
				assign=where+"disconcord";
				all.push_back("assign\t"+used_ids[id]+"\t"+assign);
				int tag=1;
				if(toi(ss1[8]) > min_frac){
					all.push_back("seq\t>"+used_ids[id]+"_read1left");
					string s=ss1[6].substr(0, toi(ss1[8]));
					all.push_back("seq\t"+s);
					all.push_back("read1left\t"+used_ids[id]+"\t"+s);
					tag=0;
				}
				if(toi(ss1[10]) > min_frac){
					all.push_back("seq\t>"+used_ids[id]+"_read1right");
					string s=ss1[6].substr(toi(ss1[8])+toi(ss1[9]), toi(ss1[10]));
					all.push_back("seq\t"+s);
					all.push_back("read1right\t"+used_ids[id]+"\t"+s);
					tag=0;
				}
				if(toi(ss2[8]) > min_frac){
					all.push_back("seq\t>"+used_ids[id]+"_read2left");
					string s=ss2[6].substr(0, toi(ss2[8]));
					all.push_back("seq\t"+s);
					all.push_back("read2left\t"+used_ids[id]+"\t"+s);
					tag=0;
				}
				if(toi(ss2[10]) > min_frac){
					all.push_back("seq\t>"+used_ids[id]+"_read2right");
					string s=ss2[6].substr(toi(ss2[8])+toi(ss2[9]), toi(ss2[10]));
					all.push_back("seq\t"+s);
					all.push_back("read2right\t"+used_ids[id]+"\t"+s);
					tag=0;
				}
				if(tag){
					all.push_back("range\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[1]+"\t"+to_string(toi(ss1[1])+toi(ss1[9]))+"\t"+ss2[0]+"\t"+ss2[1]+"\t"+to_string(toi(ss2[1])+toi(ss2[9])));
					all.push_back("assign2\t"+used_ids[id]+"\t"+assign);
					continue;
				}
				int to2=toi(ss2[1]) + toi(ss2[9]);
				int to1=toi(ss1[1]) + toi(ss1[9]);
				string xx="mark\t"+used_ids[id]+"\t"+ss1[0]+"\t"+ss1[2]+"\t"+ss1[8]+"\t"+ss1[1]+"\t"+to_string(to1)+"\t"+ss1[10]+"\t"+ss1[7]+"\t"+ss2[0]+"\t"+ss2[2]+"\t"+ss2[10]+"\t"+ss2[1]+"\t"+to_string(to2)+"\t"+ss2[8]+"\t"+ss2[7];
				all.push_back(xx);
				continue;
			}
		}
	}
}

void assigns(string dir,  string keep_chr, int min_frac, int cores=1, int strict=1, int tail=10, int block_size=500000){
	vector <string> homo_chr;
	vector <int> homo_from1;
	vector <int> homo_to1;
	vector <int> homo_from2;
	vector <int> homo_to2;

	fstream fhomo((string(dir)+"/homo.txt").c_str(),fstream::in);
	char idstr[1000];
    while(fhomo.getline(idstr,1000)){
    	vector <string> res=split(idstr,"\t");
    	homo_chr.push_back(res[0]);
    	homo_from1.push_back(toi(res[2]));
    	homo_to1.push_back(toi(res[3]));
    	homo_from2.push_back(toi(res[4]));
    	homo_to2.push_back(toi(res[5]));
    }
    fhomo.close();


    fstream fout1((string(dir)+"/temp_files/temp_mark.txt").c_str(),fstream::out);
	fstream fout2((string(dir)+"/temp_files/temp_fragment.txt").c_str(),fstream::out);
	fstream fout3((string(dir)+"/assign.txt").c_str(),fstream::out);
	fout3<<"tag\tid\tassign\n";
	fstream fout33((string(dir)+"/assign2.txt").c_str(),fstream::out);
	fout33<<"tag\tid\tassign\n";
	fstream fout4((string(dir)+"/temp_files/site_split.txt").c_str(),fstream::out);
	fout4<<"tag\tid\tchr1\tfrom\tgap\tchr2\tto\tdirection\tscore\tassign\n";
	fstream fout5((string(dir)+"/temp_files/read_split.txt").c_str(),fstream::out);
	fout5<<"tag\tid\tdirection\tchr1\tfrom1\tto1\tchr2\tfrom2\tto2\tchr3\tfrom3\tto3\tchr4\tfrom4\tto4\n";
	fstream fout6((string(dir)+"/temp_files/read_nonsplit.txt").c_str(),fstream::out);
	fout6<<"tag\tid\tchr\tdirect\tfrom\tto\n";
	fstream fout7((string(dir)+"/temp_files/temp_seq.txt").c_str(),fstream::out);
	fstream fout8((string(dir)+"/temp_files/range_split.txt").c_str(),fstream::out);
	fout8<<"tag\tid\tchr1\tfrom1\tto1\tchr2\tfrom2\tto2\n";
	fstream fout9((string(dir)+"/warnings.txt").c_str(),fstream::out);
	fstream fout10((string(dir)+"/temp_files/read_both.txt").c_str(),fstream::out);
	fout10<<"tag\tid\tchr1\tfrom\tto\n";
	fstream fout11((string(dir)+"/mapping_summary.txt").c_str(),fstream::out);
	fout11<<"tag\tNo.total\tNo.mapped\n";

	ifstream fg((string(dir)+"/aln_genome.sam.gz").c_str(), ios_base::in | ios_base::binary);
	ifstream fi((string(dir)+"/aln_insert.sam.gz").c_str(), ios_base::in | ios_base::binary);

    if(fg.good()){
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbufg;
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbufi;
		inbufg.push(boost::iostreams::gzip_decompressor());
		inbufi.push(boost::iostreams::gzip_decompressor());
		inbufg.push(fg);
		inbufi.push(fi);
		istream fgenome1(&inbufg);
		istream finsert1(&inbufi);
		while(fgenome1.getline(idstr,1000)){
		    if(idstr[0] != '@')
		    	break;
		}
		fgenome1.getline(idstr,1000);
		while(finsert1.getline(idstr,1000)){
		    if(idstr[0] != '@')
		    	break;
		}
		finsert1.getline(idstr,1000);
		int tot=0;
		char idstr1[1000];
		char idstr2[1000];
		char idstr3[1000];
		char idstr4[1000];

		char ** ary1 = new char*[block_size];
		char ** ary2 = new char*[block_size];
		char ** ary3 = new char*[block_size];
		char ** ary4 = new char*[block_size];
		for(int i = 0; i < block_size; ++i){
    		ary1[i] = new char[1000];
    		ary2[i] = new char[1000];
    		ary3[i] = new char[1000];
    		ary4[i] = new char[1000];
    	}
    	int pp=0;
    	bool efs=false;
		while(1){
		    if(pp < block_size && !efs){
		    	fgenome1.getline(idstr1,1000);
		    	fgenome1.getline(idstr2,1000);
		    	finsert1.getline(idstr3,1000);
		    	finsert1.getline(idstr4,1000);

                efs=fgenome1.eof();
		    	strcpy(ary1[pp],idstr1);
		    	strcpy(ary2[pp],idstr2);
		    	strcpy(ary3[pp],idstr3);
		    	strcpy(ary4[pp],idstr4);
		    	pp++;
		    }
		    else{
		    	if(efs)
		    		pp--;
		    	//cout<<pp<<endl;
		    	vector <string> read1;
				vector <string> read2;
				vector <string> read3;
				vector <string> read4;
				vector <string> used_ids;

                vector <string> left1;
                vector <string> left2;
		    	int nProcessors=cores;
    			omp_set_num_threads(nProcessors);
    			#pragma omp parallel for
		    	for(int j=0; j< pp; j++){
					vector <string> res1=split(ary1[j],"\t");
					vector <string> res2=split(ary2[j],"\t");
					vector <string> res3=split(ary3[j],"\t");
					vector <string> res4=split(ary4[j],"\t");
					if(res1[0] !=res2[0] || res2[0]!=res3[0] || res3[0]!=res4[0]){
					  Rcpp::Rcout<<"Error: the sam file has wrong id order!\n Please stop this program and check the 'sam' files\n";
					}
					vector <string> s1=readreads(ary1[j]);
					vector <string> s2=readreads(ary2[j]);
					vector <string> s3=readreads(ary3[j]);
					vector <string> s4=readreads(ary4[j]);

					vector <string> r1=split(s1[1]);
					vector <string> r2=split(s2[1]);
					vector <string> r3=split(s3[1]);
					vector <string> r4=split(s4[1]);
#pragma omp critical
{
					if(s1[0]=="1")
						read1.push_back(s1[1]);
					else if(s1[0]=="2")
						read2.push_back(s1[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res1[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					if(s2[0]=="1")
						read1.push_back(s2[1]);
					else if(s2[0]=="2")
						read2.push_back(s2[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res2[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					if(s3[0]=="1")
						read3.push_back(s3[1]);
					else if(s3[0]=="2")
						read4.push_back(s3[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res3[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					if(s4[0]=="1")
						read3.push_back(s4[1]);
					else if(s4[0]=="2")
						read4.push_back(s4[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res4[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					used_ids.push_back(res1[0]);
}
					r1.clear();
					r2.clear();
					r3.clear();
					r4.clear();
					s1.clear();
					s2.clear();
					s3.clear();
					s4.clear();
					res1.clear();
					res2.clear();
					res3.clear();
					res4.clear();
				}
				pp=0;

				vector <string> output;
				work(output, used_ids, read1, read2, read3, read4, homo_chr, homo_from1, homo_to1, homo_from2, homo_to2, min_frac,  strict, tail);
				int num=output.size();
				for(int i=0;i<num;i++){
					vector <string> res=split(output[i],"\t");
					if(res[0] == "mark")
						fout1<<output[i]<<endl;
					else if(res[0] == "seq")
						fout2<<output[i].substr(4)<<endl;
					else if(res[0] == "assign")
						fout3<<output[i]<<endl;
					else if(res[0] == "assign2")
						fout33<<output[i]<<endl;
					else if(res[0] == "insert")
						fout4<<output[i]<<endl;
					else if(res[0]=="split")
						fout5<<output[i]<<endl;
					else if(res[0]=="nonsplit" || res[0]=="nonsplitleft" || res[0]=="nonsplitright")
						fout6<<output[i]<<endl;
					else if(res[0] =="read1left" || res[0] =="read1right" || res[0] =="read2left" || res[0] =="read2right")
						fout7<<output[i]<<endl;
					else if(res[0]=="range")
						fout8<<output[i]<<endl;
					else if(res[0] =="warning")
						fout9<<output[i]<<endl;
					else if(res[0] =="assignboth")
						fout10<<output[i]<<endl;
					else if(res[0] =="summary")
						fout11<<output[i]<<endl;
				}
				read1.clear();
				read2.clear();
				read3.clear();
				read4.clear();
				used_ids.clear();
				output.clear();
				if(efs)
					break;
			}
		}

		for(int i = 0; i < block_size; ++i){
    		delete ary1[i];
    		delete ary2[i];
    		delete ary3[i];
    		delete ary4[i];
    	}
    	delete ary1 ;
    	delete ary2 ;
    	delete ary3 ;
    	delete ary4 ;
   	}
   	else{
		ifstream fgenome((string(dir)+"/aln_genome.sam").c_str(), ios_base::in );
		ifstream finsert((string(dir)+"/aln_insert.sam").c_str(), ios_base::in );
		streampos oldpos1 = fgenome.tellg();
		streampos oldpos2 = finsert.tellg();

		while(fgenome.getline(idstr,1000)){
		    if(idstr[0] != '@')
		    	break;
		    oldpos1 = fgenome.tellg();
		}
		fgenome.seekg(oldpos1);

		while(finsert.getline(idstr,1000)){
		    if(idstr[0] != '@')
		    	break;
		    oldpos2 = finsert.tellg();
		}
		finsert.seekg(oldpos2);

		char idstr1[1000];
		char idstr2[1000];
		char idstr3[1000];
		char idstr4[1000];
		char ** ary1 = new char*[block_size];
		char ** ary2 = new char*[block_size];
		char ** ary3 = new char*[block_size];
		char ** ary4 = new char*[block_size];
		for(int i = 0; i < block_size; ++i){
    		ary1[i] = new char[1000];
    		ary2[i] = new char[1000];
    		ary3[i] = new char[1000];
    		ary4[i] = new char[1000];
    	}
    	int pp=0;
    	bool efs=false;
		while(1){
		    if(pp < block_size && !efs){
		    	fgenome.getline(idstr1,1000);
		    	fgenome.getline(idstr2,1000);
		    	finsert.getline(idstr3,1000);
		    	finsert.getline(idstr4,1000);
		    	efs=fgenome.eof();
		    	strcpy(ary1[pp],idstr1);
		    	strcpy(ary2[pp],idstr2);
		    	strcpy(ary3[pp],idstr3);
		    	strcpy(ary4[pp],idstr4);
		    	pp++;
		    }
		    else{
		    	if(efs)
		    		pp--;
		    	vector <string> read1;
				vector <string> read2;
				vector <string> read3;
				vector <string> read4;
				vector <string> used_ids;

                vector <string> left1;
                vector <string> left2;
		    	int nProcessors=cores;
    			omp_set_num_threads(nProcessors);
    			#pragma omp parallel for
		    	for(int j=0; j< pp; j++){
					vector <string> res1=split(ary1[j],"\t");
					vector <string> res2=split(ary2[j],"\t");
					vector <string> res3=split(ary3[j],"\t");
					vector <string> res4=split(ary4[j],"\t");
					if(res1[0] !=res2[0] || res2[0]!=res3[0] || res3[0]!=res4[0]){
					  Rcpp::Rcout<<"Error: the sam file has wrong id order!\n Please stop this program and check the 'sam' files\n";
						//return(all);
					}
					vector <string> s1=readreads(ary1[j]);
					vector <string> s2=readreads(ary2[j]);
					vector <string> s3=readreads(ary3[j]);
					vector <string> s4=readreads(ary4[j]);
					vector <string> r1=split(s1[1]);
					vector <string> r2=split(s2[1]);
					vector <string> r3=split(s3[1]);
					vector <string> r4=split(s4[1]);
#pragma omp critical
{

					if(s1[0]=="1")
						read1.push_back(s1[1]);
					else if(s1[0]=="2")
						read2.push_back(s1[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res1[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					if(s2[0]=="1")
						read1.push_back(s2[1]);
					else if(s2[0]=="2")
						read2.push_back(s2[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res2[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					if(s3[0]=="1")
						read3.push_back(s3[1]);
					else if(s3[0]=="2")
						read4.push_back(s3[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res3[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					if(s4[0]=="1")
						read3.push_back(s4[1]);
					else if(s4[0]=="2")
						read4.push_back(s4[1]);
					else{
					  Rcpp::Rcout<<"Error\t"+res4[0]+": wrong reads assignment! Skipped!\n Please check 'sam' file";
					}
					used_ids.push_back(res1[0]);
}
				}
				pp=0;

				vector <string> output;
				work(output, used_ids, read1, read2, read3, read4, homo_chr, homo_from1, homo_to1, homo_from2, homo_to2, min_frac,  strict, tail);
				int num=output.size();
				for(int i=0;i<num;i++){
					vector <string> res=split(output[i],"\t");
					if(res[0] == "mark")
						fout1<<output[i]<<endl;
					else if(res[0] == "seq")
						fout2<<output[i].substr(4)<<endl;
					else if(res[0] == "assign")
						fout3<<output[i]<<endl;
					else if(res[0] == "assign2")
						fout33<<output[i]<<endl;
					else if(res[0] == "insert")
						fout4<<output[i]<<endl;
					else if(res[0]=="split")
						fout5<<output[i]<<endl;
					else if(res[0]=="nonsplit" || res[0]=="nonsplitleft" || res[0]=="nonsplitright")
						fout6<<output[i]<<endl;
					else if(res[0] =="read1left" || res[0] =="read1right" || res[0] =="read2left" || res[0] =="read2right")
						fout7<<output[i]<<endl;
					else if(res[0]=="range")
						fout8<<output[i]<<endl;
					else if(res[0] =="warning")
						fout9<<output[i]<<endl;
					else if(res[0] =="assignboth")
						fout10<<output[i]<<endl;
					else if(res[0] =="summary")
						fout11<<output[i]<<endl;
				}
				read1.clear();
				read2.clear();
				read3.clear();
				read4.clear();
				used_ids.clear();
				output.clear();
				if(efs)
					break;
			}
		}
		for(int i = 0; i < block_size; ++i){
    		delete ary1[i];
    		delete ary2[i];
    		delete ary3[i];
    		delete ary4[i];
    	}
    	delete ary1 ;
    	delete ary2 ;
    	delete ary3 ;
    	delete ary4 ;
		fgenome.close();
		finsert.close();
   	}
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
	fout6.close();
	fout7.close();
	fout8.close();
	fout9.close();
	fout10.close();
	fout11.close();
}



// [[Rcpp::export]]
void assignreads(Rcpp::StringVector ipt, Rcpp::IntegerVector para){
	// min_frac, core, strict, tail
	string dir=string(ipt(0));
	string chr=string(ipt(1));
	assigns(dir, chr, int(para(0)),int(para(1)), int(para(2)), int(para(3)), int(para(4)));
}

void readsegmentalignment(string dir, int dist=100, int block_size=500000, int unique_dif=50){
	vector <string> homo_chr;
	vector <int> homo_from1;
	vector <int> homo_to1;
	vector <int> homo_from2;
	vector <int> homo_to2;
	fstream fhomo((string(dir)+"/homo.txt").c_str(),fstream::in);
	char idstr[1000];
  while(fhomo.getline(idstr,1000)){
  	vector <string> res=split(idstr,"\t");
  	homo_chr.push_back(res[0]);
  	homo_from1.push_back(toi(res[2]));
  	homo_to1.push_back(toi(res[3]));
  	homo_from2.push_back(toi(res[4]));
  	homo_to2.push_back(toi(res[5]));
  }
  fhomo.close();

	map <string, string> gread1right, gread1left, gread2right, gread2left;
	map <string, string> iread1right, iread1left, iread2right, iread2left;
	map <string, string> seq1right, seq1left, seq2right, seq2left;
	map <string, string> mark;
	map <string, string> assign;
	map <string, int> have;
	fstream fann((string(dir)+"/temp_files/temp_seq.txt").c_str(),fstream::in);
	while(fann.getline(idstr,1000)){
        vector <string> ss=split(idstr,"\t");
        if(ss[0] =="read1left")
        	seq1left[ss[1]]=ss[2];
        else if(ss[0] =="read1right")
        	seq1right[ss[1]]=ss[2];
        else if(ss[0] =="read2left")
        	seq2left[ss[1]]=ss[2];
        else if(ss[0] =="read2right")
        	seq2right[ss[1]]=ss[2];
    }
    fann.close();

    fstream fgenome((string(dir)+"/temp_files/fragment_genome.sam").c_str(),fstream::in);
    while(fgenome.getline(idstr,1000)){
        if(idstr[0] =='@')
        	continue;
        vector <string> ss=split(idstr,"\t");
        if(ss[2] == "*")
        	continue;
        map <int,int> tg=flags(toi(ss[1]));
        if(tg[4]==1)
        	continue;
        int ll[3];
        extract_soft(ss[5],ll);
   		char *chptr;
		chptr= strchr(idstr, '_');
		int loc=chptr-idstr;
		char pid[100];
		char mytype[1000];
		memcpy( pid, &idstr[0], loc );
		pid[loc]='\0';
		string id=pid;
		int loc2=strchr(idstr, '\t')-idstr;
		memcpy( mytype, &idstr[loc+1], loc2-loc-1);
		mytype[loc2-loc-1]='\0';
		int sc=scores(idstr);
		int sc2=scores2(idstr);
		if(sc - sc2 < unique_dif)
			continue;
		string strand="0";
        have[id]=1;
		if(string(mytype) == "read1right"){
			if(ss[9] == seq1right[id])
				strand="1";
			else
				strand="-1";
			gread1right[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else if(string(mytype) == "read1left"){
			if(ss[9] == seq1left[id])
				strand="1";
			else
				strand="-1";
			gread1left[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else if(string(mytype) == "read2right"){
			if(ss[9] == seq2right[id])
				strand="1";
			else
				strand="-1";
			gread2right[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else if(string(mytype) == "read2left"){
			if(ss[9] == seq2left[id])
				strand="1";
			else
				strand="-1";
			gread2left[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else
		  Rcpp::Rcout<<"Warning: "+string(mytype)+": unknown mytype"<<endl;
	}
	fgenome.close();
	fstream finsert((string(dir)+"/temp_files/fragment_insert.sam").c_str(),fstream::in);
    while(finsert.getline(idstr,1000)){
        if(idstr[0] =='@')
        	continue;
        vector <string> ss=split(idstr,"\t");
        if(ss[2] == "*")
        	continue;

        map <int,int> tg=flags(toi(ss[1]));
        if(tg[4]==1)
        	continue;
        int ll[3];
        extract_soft(ss[5],ll);
   		//int idstr_len=strchr(idstr, '\0')-idstr;
   		char *chptr;
		chptr= strchr(idstr, '_');
		int loc=chptr-idstr;
		char pid[100];
		char mytype[1000];
		memcpy( pid, &idstr[0], loc );
		pid[loc]='\0';
		string id=pid;
		int loc2=strchr(idstr, '\t')-idstr;
		memcpy( mytype, &idstr[loc+1], loc2-loc-1);
		mytype[loc2-loc-1]='\0';
		int sc=scores(idstr);
		int sc2=scores2(idstr);
		if(sc - sc2 < 20)
			continue;
		string strand="0";
		have[id]=1;
		if(string(mytype) == "read1right"){
			if(ss[9] == seq1right[id])
				strand="1";
			else
				strand="-1";
			iread1right[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else if(string(mytype) == "read1left"){
			if(ss[9] == seq1left[id])
				strand="1";
			else
				strand="-1";
			iread1left[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else if(string(mytype) == "read2right"){
			if(ss[9] == seq2right[id])
				strand="1";
			else
				strand="-1";
			iread2right[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else if(string(mytype) == "read2left"){
			if(ss[9] == seq2left[id])
				strand="1";
			else
				strand="-1";
			iread2left[id]= ss[2]+"\t"+ss[3]+"\t"+ss[5]+"\t"+strand+"\t"+ to_string(sc)+"\t"+ss[9];
		}
		else
		  Rcpp::Rcout<<"Warning: "+string(mytype)+": unknown type"<<endl;
	}
	finsert.close();
    //cout<<"..Finish read alignment\n";

    fstream fmark((string(dir)+"/temp_files/temp_mark.txt").c_str(),fstream::in);
	vector <string> used_ids;

    map <string, int> has;
	while(fmark.getline(idstr,1000)){
        vector <string> ss=split(idstr,"\t");
        if(have[ss[1]]!=1 || has[ss[1]]==1)
            continue;
        has[ss[1]]=1;
		mark[ss[1]]=string(idstr);
		used_ids.push_back(ss[1]);
	}
	fmark.close();
	has.clear();
    //cout<<"..Finish read mark\n";

  //cout<<have.size()<<endl;
	fstream fassign((string(dir)+"/assign2.txt").c_str(),fstream::in);
    while(fassign.getline(idstr,1000)){
        vector <string> ss=split(idstr,"\t");
         if(have[ss[1]]!=1)
            continue;
        have[ss[1]]=0;
		assign[ss[1]]=ss[2];
	}
	fassign.close();
	//cout<<"..Finish read assign\n";
	//cout<<have.size()<<endl;
	//cout<<assign.size()<<endl;
	have.clear();
	ofstream fout4((string(dir)+"/temp_files/site_split.txt").c_str(),  fstream::app);
	ofstream fout5((string(dir)+"/temp_files/read_split.txt").c_str(),  fstream::app);
	ofstream fout7((string(dir)+"/temp_files/range_split.txt").c_str(), fstream::app);

	int num_used_ids=used_ids.size();
	//cout<<num_used_ids<<endl;
	int x=0;
	for(int xx=0; xx < num_used_ids; xx=xx + block_size){
		int tst=xx+ block_size;
		if(tst > num_used_ids)
			tst = num_used_ids;
		vector <string> all;
		for(int yy=xx;yy < tst; yy++){
			string id=used_ids[yy];
			//if(name != "chr17-36598226")
			//	continue;
		    //  cout<< id <<endl;

			if(iread1left[id]!="" && gread1left[id]!=""){
				vector <string> s1=split(iread1left[id],"\t");
				vector <string> s2=split(gread1left[id],"\t");
				int ll[3];
		    	extract_soft(s2[2],ll);
				int ishomo1=ishomo(s1[0], toi(s1[1]), toi(s1[1])+ ll[1], homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
				if(ishomo1)
					if(toi(s2[4]) > toi(s1[4])+20)
						iread1left[id]="";
					else
						gread1left[id]="";
				else
					if(toi(s2[4]) > toi(s1[4]))
						iread1left[id]="";
					else
						gread1left[id]="";
			}
			if(iread1right[id]!="" && gread1right[id]!=""){
				vector <string> s1=split(iread1right[id],"\t");
				vector <string> s2=split(gread1right[id],"\t");
				int ll[3];
		    	extract_soft(s2[2],ll);
				int ishomo1=ishomo(s2[0], toi(s2[1]), toi(s2[1])+ ll[1], homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
				if(ishomo1)
					if(toi(s2[4]) > toi(s1[4])+20)
						iread1right[id]="";
					else
						gread1right[id]="";
				else
					if(toi(s2[4]) > toi(s1[4]))
						iread1right[id]="";
					else
						gread1right[id]="";
			}
			if(iread2left[id]!="" && gread2left[id]!=""){
				vector <string> s1=split(iread2left[id],"\t");
				vector <string> s2=split(gread2left[id],"\t");
				int ll[3];
		    	extract_soft(s2[2],ll);
				int ishomo1=ishomo(s2[0], toi(s2[1]), toi(s2[1])+ ll[1], homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
				if(ishomo1)
					if(toi(s2[4]) > toi(s1[4])+20)
						iread2left[id]="";
					else
						gread2left[id]="";
				else
					if(toi(s2[4]) > toi(s1[4]))
						iread2left[id]="";
					else
						gread2left[id]="";
			}
			if(iread2right[id]!="" && gread2right[id]!=""){
				vector <string> s1=split(iread2right[id],"\t");
				vector <string> s2=split(gread2right[id],"\t");
				int ll[3];
		    	extract_soft(s2[2],ll);
				int ishomo1=ishomo(s2[0], toi(s2[1]), toi(s2[1])+ ll[1], homo_chr, homo_from1, homo_to1, homo_from2, homo_to2);
				if(ishomo1)
					if(toi(s2[4]) > toi(s1[4])+20)
						iread2right[id]="";
					else
						gread2right[id]="";
				else
					if(toi(s2[4]) > toi(s1[4]))
						iread2right[id]="";
					else
						gread2right[id]="";
			}
			vector <string> aa=split(mark[id],"\t");

		    //cout<<iread1left[id]<<endl;
			if(assign[id]=="split14" || assign[id]=="split23" || assign[id]=="insertdisconcord" || assign[id]=="genomedisconcord"){
				int tag=1;
				if(iread1left[id]!="" && assign[id]!="split23" && assign[id]!="genomedisconcord"){
					vector <string> ss=split(iread1left[id],"\t");
					//cout<<id<<" "<<iread1left[id]<<endl;
					//cout<<id<<" "<<mark[id]<<endl;
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq1left[id])
						ff="rf";
					if(ss[0] != aa[9])
						x=1;
					else{
						int gap=toi(ss[1])-toi(aa[13]);
						int from=toi(ss[1])+ll[1];
						int dd=ll[2];
						if(ff=="fr" || ff=="rf"){
							gap=toi(aa[12])- toi(ss[1]) - ll[1];
							from=toi(ss[1]);
							dd=ll[0];
						}
						//cout<<ff<<" "<<gap<<endl;
						if(gap < dist && gap > -20){
							bool is=msc(ss[0],to_string(from),aa[2],aa[5],ff);
							all.push_back(ms(is, id,ss[0],to_string(from),to_string(dd),aa[2],aa[5],ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))), assign[id]+"1left"));
							all.push_back(mc(is, id,ff,aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],"","",""));
							//cout<<"1: "<<mc(id,ff,aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],"","","")<<endl;
							tag=0;
						}
					}
				}
				if(iread1right[id]!="" && assign[id]!="split23" && assign[id]!="genomedisconcord"){
					vector <string> ss=split(iread1right[id],"\t");
					int ll[3];
		   			extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq1right[id])
						ff="fr";
					if(ss[0] != aa[9])
						x=1;
					else{
						int gap=toi(aa[12])- toi(ss[1]) - ll[1];
						int to=toi(ss[1]);
						int dd=ll[0];
						if(ff=="fr" || ff=="rf"){
							gap=toi(ss[1])-toi(aa[13]);
							to=toi(ss[1])+ll[1];
							dd=ll[2];
						}
						if(gap < dist && gap > -20){
							bool is=msc(aa[2],aa[6],ss[0],to_string(to),ff);
							all.push_back(ms(is, id,aa[2],aa[6],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))),assign[id]+"1right"));
							all.push_back(mc(is, id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[7],aa[12],aa[13]));
							//cout<<"2: "<<mc(id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[7],aa[12],aa[13])<<endl;
							tag=0;
						}
					}
				}
				if(iread2left[id]!="" && assign[id]!="split14" && assign[id]!="genomedisconcord"){
					vector <string> ss=split(iread2left[id],"\t");
					int ll[3];
		   			extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq2left[id])
						ff="rf";
					if(ss[0] != aa[2])
						x=1;
					else{
						int gap=toi(ss[1])-toi(aa[6]);
						int dd=ll[2];
						int from=toi(ss[1])+ll[1];
						if(ff=="fr"|| ff=="rf"){
							gap=toi(aa[5])- toi(ss[1]) - ll[1];
							dd=ll[0];
							from=toi(ss[1]);
						}
						if(gap < dist && gap > -20){
							bool is=msc(ss[0],to_string(from),aa[9],aa[12],ff);
							all.push_back(ms(is,id,ss[0],to_string(from),to_string(dd),aa[9],aa[12],ff, to_string(sqrt(toi(ss[4]) * toi(aa[15]))),assign[id]+"2left"));
							all.push_back(mc(is,id,ff,aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[9],aa[12],aa[13],"","",""));
							//cout<<"3: "<<mc(id,ff,"","","",aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6])<<endl;
							tag=0;
						}
					}
				}
				if(iread2right[id]!="" && assign[id]!="split14" && assign[id]!="genomedisconcord"){
					vector <string> ss=split(iread2right[id],"\t");
					int ll[3];
		   			 extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq2right[id])
						ff="fr";
					if(ss[0] != aa[2])
						x=1;
					else{
						int gap=toi(aa[5])- toi(ss[1]) - ll[1];
						int dd=ll[0];
						int to=toi(ss[1]);
						if(ff=="fr"|| ff=="rf"){
							gap=toi(ss[1])-toi(aa[6]);
							dd=ll[2];
							to=toi(ss[1])+ll[1];
						}
						if(gap < dist && gap > -20){
							bool is =msc(aa[9],aa[13],ss[0],to_string(to),ff);
							all.push_back(ms(is, id,aa[9],aa[13],to_string(dd),ss[0],to_string(to),ff, to_string(sqrt(toi(ss[4]) * toi(aa[15]))),assign[id]+"2right"));
							all.push_back(mc(is, id,ff,"","","",aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6]));
							//cout<<"4: "<<mc(id,ff,"","","",aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6])<<endl;
							tag=0;
						}
					}
				}
				if(gread1left[id]!="" && assign[id]!="split14" && assign[id]!="insertdisconcord"){
					vector <string> ss=split(gread1left[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq1left[id])
						ff="rf";
					if(ss[0] != aa[9])
						x=1;
					else{
						int gap=toi(ss[1])-toi(aa[13]);
						int from=toi(ss[1])+ll[1];
						int dd=ll[2];
						if(ff=="fr"|| ff=="rf"){
							gap=toi(aa[12])- toi(ss[1]) - ll[1];
							from=toi(ss[1]);
							dd=ll[0];
						}
						if(gap < dist && gap > -20){
							bool is =msc(ss[0],to_string(from),aa[2],aa[5],ff);
							all.push_back(ms(is,id,ss[0],to_string(from),to_string(dd),aa[2],aa[5],ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))),assign[id]+"1left"));
							all.push_back(mc(is,id,ff,aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],"","",""));
							//cout<<"5: "<<mc(id,ff,aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],"","","")<<endl;
							tag=0;
						}
					}
				}
				if(gread1right[id]!="" && assign[id]!="split14" && assign[id]!="insertdisconcord"){
					vector <string> ss=split(gread1right[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq1right[id])
						ff="fr";
					if(ss[0] != aa[9])
						x=1;
					else{
						int gap=toi(aa[12])- toi(ss[1]) - ll[1];
						int to=toi(ss[1]);
						int dd=ll[0];
						if(ff=="fr"|| ff=="rf"){
							gap=toi(ss[1])-toi(aa[13]);
							to=toi(ss[1])+ll[1];
							dd=ll[2];
						}
						if(gap < dist && gap > -20){
							bool is=msc(aa[2],aa[6],ss[0],to_string(to),ff);
							all.push_back(ms(is,id,aa[2],aa[6],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))),assign[id]+"1right"));
							all.push_back(mc(is,id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[9],aa[12],aa[13]));
							//cout<<"6: "<<mc(id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[7],aa[12],aa[13])<<endl;
							tag=0;
						}
					}
				}
				if(gread2left[id]!="" && assign[id]!="split23" && assign[id]!="insertdisconcord"){
					vector <string> ss=split(gread2left[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="rf";
					if(ss[5] != seq2left[id])
						ff="fr";
					if(ss[0] != aa[2])
						x=1;
					else{
						int gap=toi(ss[1])-toi(aa[6]);
						int from=toi(ss[1])+ll[1];
						int dd=ll[2];
						if(ff=="fr"|| ff=="rf"){
							gap=toi(aa[5])- toi(ss[1]) - ll[1];
							from=toi(ss[1]);
							dd=ll[0];
						}
						if(gap < dist && gap > -20){
							bool is=msc(ss[0],to_string(from),aa[9],aa[12],ff);
							all.push_back(ms(is, id,ss[0],to_string(from),to_string(dd),aa[9],aa[12],ff, to_string(sqrt(toi(ss[4]) * toi(aa[15]))),assign[id]+"2left"));
							all.push_back(mc(is, id,ff,"","","",aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6]));
							//cout<<"7: "<<mc(id,ff,"","","",aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6])<<endl;
							tag=0;
						}
					}
				}
				if(gread2right[id]!="" && assign[id]!="split23" && assign[id]!="insertdisconcord"){
					vector <string> ss=split(gread2right[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq2right[id])
						ff="fr";
					if(ss[0] != aa[2])
						x=1;
					else{
						int gap=toi(aa[5])- toi(ss[1]) - ll[1];
						int to=toi(ss[1]);
						int dd=ll[0];
						if(ff=="fr"|| ff=="rf"){
							gap=toi(ss[1])-toi(aa[6]);
							to=toi(ss[1])+ll[1];
							dd=ll[2];
						}
						if(gap < dist && gap > -20){
							bool is=msc(aa[9],aa[13],ss[0],to_string(to),ff);
							all.push_back(ms(is,id,aa[9],aa[13],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[15]))),assign[id] +"2right"));
							all.push_back(mc(is,id,ff,"","","",aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6]));
							//cout<<"7: "<<mc(id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[7],aa[12],aa[13])<<endl;
							tag=0;
						}
					}
				}
				if(tag){
					all.push_back("range\t"+id+"\t"+aa[2]+"\t"+aa[5]+"\t"+aa[6]+"\t"+aa[9]+"\t"+aa[12]+"\t"+aa[13]);
				}
			}
			else if(assign[id]=="insertconcord" || assign[id]=="genomeconcord"){
				int tag=1;
				if(iread1left[id]!=""){
					vector <string> ss=split(iread1left[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq1left[id])
						ff="rf";
					int from=toi(ss[1])+ll[1];
					int dd=ll[2];
					if(ff=="fr" || ff=="rf"){
						from=toi(ss[1]);
						dd=ll[0];
					}
					int gap=10000;
					if(ss[0]==aa[2])
						gap=abs(from-toi(aa[5]));
					//cout<<ff<<" "<< dd << " "<< gap <<" " << ll[0]<<" " << ll[2]<<" "<<ss[0]<<"/"<<aa[2] <<endl;
					if(!(ll[0] > 20 && ll[2] > 20) ){
						bool is=msc(ss[0],to_string(from),aa[2],aa[5],ff);
						all.push_back(ms(is,id,ss[0],to_string(from),to_string(dd),aa[2],aa[5],ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))),assign[id]+"1left"));
						if(aa.size()>9)
							all.push_back(mc(is,id,ff,"","","",ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],aa[9],aa[12],aa[13]));
						else
							all.push_back(mc(is,id,ff,"","","",ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],"","",""));
						tag=0;
					}
				}
				if(iread2right[id]!=""){
					vector <string> ss=split(iread2right[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq2right[id])
						ff="fr";

					int to=toi(ss[1]);
					int dd=ll[0];
					if(ff=="fr"|| ff=="rf"){
						to=toi(ss[1])+ll[1];
						dd=ll[2];
					}
					int gap=10000;
					if(aa[2]==ss[0])
						gap=abs(to-toi(aa[6]));
					if(aa.size() > 9){
						if(!(ll[0] > 20 && ll[2] > 20) ){
							bool is=msc(aa[9],aa[13],ss[0],to_string(to),ff);
							all.push_back(ms(is,id,aa[9],aa[13],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[15]))),assign[id]+"2right"));
							all.push_back(mc(is,id, ff,aa[2],aa[5],aa[6],aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),"","",""));
							tag=0;
						}
					}
					else{
						if(!(ll[0] > 20 && ll[2] > 20) ){
							bool is=msc(aa[2],aa[6],ss[0],to_string(to),ff);
							all.push_back(ms(is,id,aa[2],aa[6],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))),assign[id]+"2right"));
							all.push_back(mc(is,id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),"","",""));
							tag=0;
						}
					}
				}
				if(gread1left[id]!=""){
					vector <string> ss=split(gread1left[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq1left[id])
						ff="rf";
					int from=toi(ss[1])+ll[1];
					int dd=ll[2];
					if(ff=="fr" || ff=="rf"){
						from=toi(ss[1]);
						dd=ll[0];
					}
					int gap=10000;
					if(ss[0]==aa[2])
						gap=abs(from-toi(aa[5]));
					if(!(ll[0] > 20 && ll[2] > 20)){
						bool is=msc(ss[0],to_string(from),aa[2],aa[5],ff);
						all.push_back(ms(is,id,ss[0],to_string(from),to_string(dd),aa[2],aa[5],ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))),assign[id]+"1left"));
						if(aa.size()>9)
							all.push_back(mc(is,id,ff,"","","",ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],aa[9],aa[12],aa[13]));
						else
							all.push_back(mc(is,id,ff,"","","",ss[0],ss[1],to_string(toi(ss[1])+ll[1]),aa[2],aa[5],aa[6],"","",""));
						tag=0;
					}
				}
				if(gread2right[id]!=""){
					vector <string> ss=split(gread2right[id],"\t");
					int ll[3];
		    		extract_soft(ss[2],ll);
					string ff="ff";
					if(ss[5] != seq2right[id])
						ff="fr";
					int to=toi(ss[1]);
					int dd=ll[0];
					if(ff=="fr"|| ff=="rf"){
						to=toi(ss[1])+ll[1];
						dd=ll[2];
					}
					int gap=10000;
					if(aa[2]==ss[0])
						gap=abs(to-toi(aa[6]));
					if(aa.size() > 9){
						if(!( ll[0] > 20 && ll[2] > 20) ){
							bool is=msc(aa[9],aa[13],ss[0],to_string(to),ff);
							all.push_back(ms(is,id,aa[9],aa[13],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[15]))),assign[id]+"2right"));
							all.push_back(mc(is,id, ff,aa[2],aa[5],aa[6],aa[9],aa[12],aa[13],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),"","",""));
							tag=0;
						}
					}
					else{
						if(!( ll[0] > 20 && ll[2] > 20 ) ){
							bool is=msc(aa[2],aa[6],ss[0],to_string(to),ff);
							all.push_back(ms(is,id,aa[2],aa[6],to_string(dd),ss[0],to_string(to),ff,to_string(sqrt(toi(ss[4]) * toi(aa[8]))), assign[id]+"2right"));
							all.push_back(mc(is,id,ff,"","","",aa[2],aa[5],aa[6],ss[0],ss[1],to_string(toi(ss[1])+ll[1]),"","",""));
							tag=0;
						}
					}
				}
				if(tag){
					if(aa.size() > 9)
						all.push_back("range\t"+id+"\t"+aa[2]+"\t"+aa[5]+"\t"+aa[6]+"\t"+aa[9]+"\t"+aa[12]+"\t"+aa[13]);
					else
						all.push_back("range\t"+id+"\t"+aa[2]+"\t"+aa[5]+"\t"+aa[6]+"\tempty\t0\t0");
				}
			}
		}
		int num=all.size();
		for(int i=0;i<num;i++){
			vector <string> res=split(all[i],"\t");
			if(res[0] == "insert")
				fout4<<all[i]<<endl;
			else if(res[0]=="split")
				fout5<<all[i]<<endl;
			else if(res[0]=="range")
				fout7<<all[i]<<endl;
		}
	}
	fout4.close();
	fout5.close();
	fout7.close();
}


// [[Rcpp::export]]
void reassignreads(Rcpp::String dir, Rcpp::IntegerVector dist){
	readsegmentalignment(dir, int(dist(0)), int(dist(1)), int(dist(2)));
}



