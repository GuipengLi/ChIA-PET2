// ChIA-PET2
// Copyright (c) 2016 Guipeng Li                               
// Author(s): Guipeng Li
// Contact: guipeng.lee@gmail.com


#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <sstream>
#include <bitset>
#include <map>
#include <iterator>
#include <stdlib.h>
#include <climits>
#include <iomanip>
//#include "full_matrix_aligner.h"
#include <zlib.h>
//#include "kseq.h"
#include <omp.h>
//KSEQ_INIT(gzFile, gzread)

using namespace std;
static const size_t npos = -1;

size_t myfind(string data, string needle){
	const char *result = strstr(data.c_str(), needle.c_str());
	size_t found;
    if (result == NULL) {
		found = -1;
	} else {																						                             found = result - data.c_str();
	}                  
	return (found);
}

string vector_join( const vector<string>& v, const string& token ){
    ostringstream result;
    for (vector<string>::const_iterator i = v.begin(); i != v.end(); i++){
        if (i != v.begin()) result << token;
        result << *i;
    }
    return result.str();
}


vector<string> string_split( const string& s, const string& delimiter ){
    vector<string> result;
    string::size_type from = 0;
    string::size_type to = 0;
    while ( to != string::npos ){
        to = s.find( delimiter, from );
        if ( from < s.size() && from != to ){
            result.push_back( s.substr( from, to - from ) );
        }
        from = to + delimiter.size();
    }
    return result;
}

string get_strand( unsigned long x ) {
    string strand = "+";
    if ( x & 0x10 )
    {
        strand = "-";
    }
    return strand;
}


int StringToInt( string Text ) {
    int output;
    if ( ! (istringstream(Text) >> output) ) output = 0;
    return output;
}

string IntToString( int Number ) {
    string Result;          // string which will contain the result
    ostringstream convert;   // stream used for the conversion
    convert << Number;      // insert the textual representation of 'Number' in the characters in the stream
    Result = convert.str(); // set 'Result' to the contents of the stream
    return Result;
}

template <typename T>
string NumberToString ( T Number )
{
  stringstream ss;
	ss << Number;
	return ss.str();
}

size_t bitap_search(string text, string pattern, int k)
{
     size_t result = npos;
     int m = pattern.length();
     unsigned long *R;
     unsigned long pattern_mask[CHAR_MAX+1];
     int i, d;
 
     if (pattern[0] == '\0') return -1;
     if (m > 31) cout << "The pattern is too long!";
 
     /* Initialize the bit array R */
     R = (unsigned long *) malloc((k+1) * sizeof(*R));
     for (i=0; i <= k; ++i)
         R[i] = ~1;
 
     /* Initialize the pattern bitmasks */
     for (i=0; i <= CHAR_MAX; ++i)
         pattern_mask[i] = ~0;
     for (i=0; i < m; ++i)
         pattern_mask[pattern[i]] &= ~(1UL << i);
 
     for (i=0; text[i] != '\0'; ++i) {
         /* Update the bit arrays */
         unsigned long old_Rd1 = R[0];
 
         R[0] |= pattern_mask[text[i]];
         R[0] <<= 1;
 
         for (d=1; d <= k; ++d) {
             unsigned long tmp = R[d];
             /* Substitution is all we care about */
             R[d] = (old_Rd1 & (R[d] | pattern_mask[text[i]])) << 1;
             old_Rd1 = tmp;
         }
 
         if (0 == (R[k] & (1UL << m))) {
             result = (i - m) + 1;
             break;
         }
     }
 
     free(R);
     return result;
}

//int findlinker(const char* target,const char* query, float error) {
//  int qlen = strlen(query);
//  int tlen = strlen(target);
//  charvec path;
//  kv_init(path);
//  result res = align_full_matrix( target, query, qlen, tlen, &path, 0);
//  //printf("Target: %d - %d\n", res.tstart, res.tend);
//  int pos = res.tstart;
//  int path_size = kv_size(path);
//  unsigned char cigar[path_size];
//  int i;
//  for(i = 0; i < path_size; i++) {
//    cigar[i] = kv_pop(path);
//  }
//  float ERR = 1 - cigar_accuracy(cigar, path_size);
//  kv_destroy(path);
//  if (ERR>error) pos = -1;
//  return(pos);
//}

/*************
 * FASTQ I/O *
 *************/

typedef struct {
	string name, seq, qual;
} bseq1_t;


vector <bseq1_t> read_fastq(gzFile &file, int chunk_size)
{
	int n=0;
	bseq1_t s;
	vector <bseq1_t> seqs;
	int LEN = 1000;
    char *buffer = new char[LEN];
	while (0!=gzgets(file, buffer, LEN)) {
		s.name = strtok(buffer,"\n");
		gzgets(file, buffer, LEN);
		if (strlen(buffer)>1)
			s.seq = strtok(buffer,"\n");
		else
			s.seq = "";
		gzgets(file, buffer, LEN);
		gzgets(file, buffer, LEN);
		if (strlen(buffer)>1)
			s.qual = strtok(buffer,"\n");
		else
			s.qual = "";
		seqs.push_back(s);
		n++;
		if (n >= chunk_size) break;
	}
	return seqs;
}



int parseFastqgz_p(int threadN,string fn1, string fn2,string basename, size_t minlength = 15,
              int keepempty = 0, string linker1 = "GTTGGATAAG" , string linker2 = "GTTGGAATGT", int msk=0)
{
	omp_set_num_threads(threadN);
    ofstream same1 ( (basename + "_1.valid.fastq").c_str() );
    ofstream same2 ( (basename + "_2.valid.fastq").c_str() );
    ofstream trimstat ( (basename + ".trim.stat").c_str() );
	gzFile fp1, fp2;
	fp1 = gzopen(fn1.c_str(), "rb");
	fp2 = gzopen(fn2.c_str(), "rb");
	int n_seqs1 = 1;
	long long linecount = 0;
	int expectcount = 0;
	int validcount = 0;
	int chimcount = 0;
	int ambicount = 0;
	int ebothcount = 0;
	int emptycount1 = 0;
	int emptycount2 = 0;
	int chunk_line = 2000000;
	vector<bseq1_t> seqs1;
	vector<bseq1_t> seqs2;
	while(n_seqs1>0){
		seqs1 = read_fastq(fp1, chunk_line);
		seqs2 = read_fastq(fp2, chunk_line);
		n_seqs1 = seqs1.size();
		#pragma omp parallel for reduction(+:linecount) reduction(+:expectcount) reduction(+:ebothcount) reduction(+:validcount) reduction(+:chimcount) reduction(+:emptycount1) reduction(+:emptycount2)
		for (int kk=0;kk<n_seqs1;kk++){
			linecount++;
			bseq1_t s1 = seqs1[kk];
			bseq1_t s2 = seqs2[kk];
			string seq1 = s1.seq;
			string qual1 = s1.qual;
			string name1 = s1.name;
			string seq2 = s2.seq;
			string qual2 = s2.qual;
			string name2 = s2.name;
			size_t s1l1=npos, s1l2=npos, s2l1=npos, s2l2=npos;
			//s1l1 = seq1.find(linker1);
			//s1l2 = seq1.find(linker2);
			//s2l1 = seq2.find(linker1);
			//s2l2 = seq2.find(linker2);
			s1l1 = myfind(seq1,linker1);
			s1l2 = myfind(seq1,linker2);
			s2l1 = myfind(seq2,linker1);
			s2l2 = myfind(seq2,linker2);
			if (msk>0){
 				if (s1l1==npos ) s1l1 = bitap_search(seq1, linker1, msk);	
 				if (s1l2==npos ) s1l2 = bitap_search(seq1, linker2, msk);	
 				if (s2l1==npos ) s2l1 = bitap_search(seq2, linker1, msk);	
 				if (s2l2==npos ) s2l2 = bitap_search(seq2, linker2, msk);	
			}
			//cout << "line "<< linecount++<<endl;
			//if (s1l1==0 || s1l2==0 || s2l1==0 || s2l2==0)
			//{cout<<"zero..."<<endl; continue;}
			int link1 = 0;		
			if (s1l1==npos && s1l2==npos)
				link1 = 0;
			else if ( s1l1!=npos && s1l2==npos){
				link1 = 1;
				seq1 = seq1.substr(0,s1l1);
				qual1 = qual1.substr(0,s1l1);
			}else if ( s1l1==npos && s1l2!=npos){
				link1 = 2;
				seq1 = seq1.substr(0,s1l2);
				qual1 = qual1.substr(0,s1l2);
			} else link1 = 3;

			int link2 = 0;		
			if (s2l1==npos && s2l2==npos)
				link2 = 0;
			else if ( s2l1!=npos && s2l2==npos){
				link2 = 1;
				seq2 = seq2.substr(0,s2l1);
				qual2 = qual2.substr(0,s2l1);
			}else if ( s2l1==npos && s2l2!=npos){
				link2 = 2;
				seq2 = seq2.substr(0,s2l2);
				qual2 = qual2.substr(0,s2l2);
			} else link2 = 3;


			//determine pairtype
			string pairtype = "ambi";
			if ( (link1==1 && link2==1) || (link1==2 && link2==2) ) {
				pairtype = "expect";
				ebothcount++;
			} else if ( (link1==1 && link2==2) || (link1==2 && link2==1) )
				pairtype = "chim";
		
			if ( keepempty==2){
				if ((link1 == 0 && link2 == 1) ||
            		(link1 == 0 && link2 == 2) ||
                	(link1 == 1 && link2 == 0) ||
                	(link1 == 2 && link2 == 0) ){
					pairtype = "expect";
					emptycount1++;
				}
            	if (link1 == 0 && link2 == 0) {
					pairtype = "expect";
					emptycount2++;
				}
			} else if (keepempty == 1) {
        		if ((link1 == 0 && link2 == 1) ||
            		(link1 == 0 && link2 == 2) ||
                	(link1 == 1 && link2 == 0) ||
                	(link1 == 2 && link2 == 0) ) {
                	pairtype = "expect";
 		    		emptycount1++;
            	}
            	if (link1 == 0 && link2 == 0) {
					emptycount2++;
				}
			} else {
        		if ((link1 == 0 && link2 == 1) ||
            		(link1 == 0 && link2 == 2) ||
                	(link1 == 1 && link2 == 0) ||
                	(link1 == 2 && link2 == 0) ) {
 		    		emptycount1++;
            	}
            	if (link1 == 0 && link2 == 0) {
					emptycount2++;
				}
			}

			if (pairtype == "expect")
				expectcount++;
        	else if (pairtype == "chim")
				chimcount++;
			else 
				ambicount++;

			#pragma omp critical
			if ( (pairtype=="expect") &&  ( seq1.length() >= minlength ) && (seq2.length() >= minlength ))
			{
				validcount++;
				same1 << name1<<endl << seq1<<endl << "+" <<endl<< qual1<<endl;
				same2 << name2<<endl << seq2<<endl << "+" <<endl<< qual2<<endl;
			}

		}
		if(n_seqs1>0) cout << "Processed "<< n_seqs1<< " pair reads"<<endl<< flush;
		seqs1.clear();
		seqs2.clear();
	}
	gzclose(fp1);
	gzclose(fp2);
	same1.close();
	same2.close();
    cout << "Total PETs:\t" << linecount <<endl;
    cout << "Expect PETs:\t" << expectcount <<endl;
    cout << "Expect both PETs:\t" << ebothcount <<endl;
    cout << "Chim PETs:\t" << chimcount <<endl;
    cout << "1Empty PETs:\t" << emptycount1 <<endl;
    cout << "2Empty PETs:\t" << emptycount2 <<endl;
    cout << "Valid PETs:\t" << validcount <<endl <<flush;
    trimstat << "Total PETs:\t" << linecount <<endl;
    trimstat << "Expect PETs:\t" << expectcount <<endl;
    trimstat << "Expect both PETs:\t" << ebothcount <<endl;
    trimstat << "Chim PETs:\t" << chimcount <<endl;
    trimstat << "1Empty PETs:\t" << emptycount1 <<endl;
    trimstat << "2Empty PETs:\t" << emptycount2 <<endl;
    trimstat << "Valid PETs:\t" << validcount <<endl <<flush;
	trimstat.close();
	return 0;
}


int parseFastqBridgegz_p(int threadN,string fn1, string fn2,string basename, size_t minlength = 15,
              int keepempty = 0, string linker1="ACGCGATATCTTATCTGACT", string linker2="AGTCAGATAAGATATCGCGT", int msk=0)
{
	omp_set_num_threads(threadN);
    ofstream same1 ( (basename + "_1.valid.fastq").c_str() );
    ofstream same2 ( (basename + "_2.valid.fastq").c_str() );
    ofstream trimstat ( (basename + ".trim.stat").c_str() );
	gzFile fp1, fp2;
	fp1 = gzopen(fn1.c_str(), "rb");
	fp2 = gzopen(fn2.c_str(), "rb");
	int n_seqs1 = 1;
	long long linecount = 0;
	int expectcount = 0;
	int validcount = 0;
	int ebothcount = 0;
	int chimcount = 0;
	int ambicount = 0;
	int emptycount1 = 0;
	int emptycount2 = 0;
	int chunk_line = 2000000;
	vector<bseq1_t> seqs1;
	vector<bseq1_t> seqs2;
	while(n_seqs1>0){
		seqs1 = read_fastq(fp1, chunk_line);
		seqs2 = read_fastq(fp2, chunk_line);
		n_seqs1 = seqs1.size();
		#pragma omp parallel for reduction(+:linecount) reduction(+:expectcount) reduction(+:ebothcount) reduction(+:validcount) reduction(+:chimcount) reduction(+:emptycount1) reduction(+:emptycount2)
		for (int kk=0;kk<n_seqs1;kk++){
			linecount++;
			bseq1_t s1 = seqs1[kk];
			bseq1_t s2 = seqs2[kk];
			string seq1 = s1.seq;
			string qual1 = s1.qual;
			string name1 = s1.name;
			string seq2 = s2.seq;
			string qual2 = s2.qual;
			string name2 = s2.name;
			size_t s1l1=npos, s1l2=npos, s2l1=npos, s2l2=npos;
			//s1l1 = seq1.find(linker1);
			//s1l2 = seq1.find(linker2);
			//s2l1 = seq2.find(linker1);
			//s2l2 = seq2.find(linker2);
			s1l1 = myfind(seq1,linker1);
			s1l2 = myfind(seq1,linker2);
			s2l1 = myfind(seq2,linker1);
			s2l2 = myfind(seq2,linker2);
			if (msk>0){
 				if (s1l1==npos ) s1l1 = bitap_search(seq1, linker1, msk);	
 				if (s1l2==npos ) s1l2 = bitap_search(seq1, linker2, msk);	
 				if (s2l1==npos ) s2l1 = bitap_search(seq2, linker1, msk);	
 				if (s2l2==npos ) s2l2 = bitap_search(seq2, linker2, msk);	
			}
			int link1 = 0;		
			if (s1l1==npos && s1l2==npos)
				link1 = 0;
			else if ( s1l1!=npos && s1l2==npos){
				link1 = 1;
				seq1 = seq1.substr(0,s1l1);
				qual1 = qual1.substr(0,s1l1);
			}else if ( s1l1==npos && s1l2!=npos){
				link1 = 2;
				seq1 = seq1.substr(0,s1l2);
				qual1 = qual1.substr(0,s1l2);
			} else link1 = 3;

			int link2 = 0;		
			if (s2l1==npos && s2l2==npos)
				link2 = 0;
			else if ( s2l1!=npos && s2l2==npos){
				link2 = 1;
				seq2 = seq2.substr(0,s2l1);
				qual2 = qual2.substr(0,s2l1);
			}else if ( s2l1==npos && s2l2!=npos){
				link2 = 2;
				seq2 = seq2.substr(0,s2l2);
				qual2 = qual2.substr(0,s2l2);
			} else link2 = 3;


			//determine pairtype
			string pairtype = "ambi";
			if ( (link1==1 && link2==2) || (link1==2 && link2==1) ){
				pairtype = "expect";
				ebothcount++;
			} else if ( (link1==1 && link2==1) || (link1==2 && link2==2) )
				pairtype = "chim";
		
			if ( keepempty==2){
				if ((link1 == 0 && link2 == 1) ||
            		(link1 == 0 && link2 == 2) ||
                	(link1 == 1 && link2 == 0) ||
                	(link1 == 2 && link2 == 0) ) {
					pairtype = "expect";
					emptycount1++;
				}
            	if (link1 == 0 && link2 == 0) {
					pairtype = "expect";
					emptycount2++;
				}
			}
			else if (keepempty == 1) {
        		if ((link1 == 0 && link2 == 1) ||
            		(link1 == 0 && link2 == 2) ||
                	(link1 == 1 && link2 == 0) ||
                	(link1 == 2 && link2 == 0) ) {
                	pairtype = "expect";
 		    		emptycount1++;
            	}
            	if (link1 == 0 && link2 == 0) {
					emptycount2++;
				}
			} else {
        		if ((link1 == 0 && link2 == 1) ||
            		(link1 == 0 && link2 == 2) ||
                	(link1 == 1 && link2 == 0) ||
                	(link1 == 2 && link2 == 0) ) {
 		    		emptycount1++;
            	}
            	if (link1 == 0 && link2 == 0) {
					emptycount2++;
				}
			}
		
			if (pairtype == "expect")
				expectcount++;
        	else if (pairtype == "chim")
				chimcount++;
			else 
				ambicount++;

			#pragma omp critical
			if ( (pairtype=="expect") &&  ( seq1.length() >= minlength ) && (seq2.length() >= minlength ))
			{
				validcount++;
				same1 << name1<<endl << seq1<<endl << "+" <<endl<< qual1<<endl;
				same2 << name2<<endl << seq2<<endl << "+" <<endl<< qual2<<endl;
			}

		}
		if (n_seqs1>0) cout << "Processed "<< n_seqs1<< " pair reads"<<endl<<flush;
		seqs1.clear();
		seqs2.clear();
	}
	gzclose(fp1);
	gzclose(fp2);
	same1.close();
	same2.close();
    cout << "Total PETs:\t" << linecount <<endl;
    cout << "Expect PETs:\t" << expectcount <<endl;
    cout << "Expect both PETs:\t" << ebothcount <<endl;
    cout << "Chim PETs:\t" << chimcount <<endl;
    //cout << "Ambi PETs: " << ambicount <<endl;
    cout << "1Empty PETs:\t" << emptycount1 <<endl;
    cout << "2Empty PETs:\t" << emptycount2 <<endl;
    cout << "Valid PETs:\t" << validcount <<endl<<flush;
    trimstat << "Total PETs:\t" << linecount <<endl;
    trimstat << "Expect PETs:\t" << expectcount <<endl;
    trimstat << "Expect both PETs:\t" << ebothcount <<endl;
    trimstat << "Chim PETs:\t" << chimcount <<endl;
    trimstat << "1Empty PETs:\t" << emptycount1 <<endl;
    trimstat << "2Empty PETs:\t" << emptycount2 <<endl;
    trimstat << "Valid PETs:\t" << validcount <<endl <<flush;
	trimstat.close();
	return 0;
}




int parseFastqEnzymegz_p(int threadN,string fn1, string fn2,string basename, size_t minlength = 15,
              int keepempty = 0, string linker1 = "GATCGATC" , string linker2 = "GATCGATC", int msk=0)
{
	omp_set_num_threads(threadN);
    ofstream same1 ( (basename + "_1.valid.fastq").c_str() );
    ofstream same2 ( (basename + "_2.valid.fastq").c_str() );
    ofstream trimstat ( (basename + ".trim.stat").c_str() );
	gzFile fp1, fp2;
	fp1 = gzopen(fn1.c_str(), "rb");
	fp2 = gzopen(fn2.c_str(), "rb");
	int n_seqs1 = 1;
	long long linecount = 0;
	int expectcount = 0;
	int ebothcount = 0;
	int validcount = 0;
	int chimcount = 0;
	int ambicount = 0;
	int emptycount1 = 0;
	int emptycount2 = 0;
	int chunk_line = 2000000;
	vector<bseq1_t> seqs1;
	vector<bseq1_t> seqs2;
	while(n_seqs1>0){
		seqs1 = read_fastq(fp1, chunk_line);
		seqs2 = read_fastq(fp2, chunk_line);
		n_seqs1 = seqs1.size();
		#pragma omp parallel for reduction(+:linecount) reduction(+:expectcount) reduction(+:ebothcount) reduction(+:validcount) reduction(+:chimcount) reduction(+:emptycount1) reduction(+:emptycount2)
		for (int kk=0;kk<n_seqs1;kk++){
			linecount++;
			bseq1_t s1 = seqs1[kk];
			bseq1_t s2 = seqs2[kk];
			string seq1 = s1.seq;
			string qual1 = s1.qual;
			string name1 = s1.name;
			string seq2 = s2.seq;
			string qual2 = s2.qual;
			string name2 = s2.name;
			size_t s1l1=npos, s2l1=npos;
			//s1l1 = seq1.find(linker1);
			//s2l1 = seq2.find(linker1);
			s1l1 = myfind(seq1,linker1);
			s2l1 = myfind(seq2,linker1);
 			if (s1l1==npos && msk>0) s1l1 = bitap_search(seq1, linker1, msk);	
 			if (s2l1==npos && msk>0) s2l1 = bitap_search(seq2, linker1, msk);	
			//cout << "line "<< linecount++<<endl;
			//if (s1l1==0 || s1l2==0 || s2l1==0 || s2l2==0)
			//{cout<<"zero..."<<endl; continue;}
			int link1 = 0;		
			if (s1l1!=npos) {
				link1 = 1;
				seq1 = seq1.substr(0,s1l1);
				qual1 = qual1.substr(0,s1l1);
			}

			int link2 = 0;		
			if (s2l1!=npos) {
				link2 = 1;
				seq2 = seq2.substr(0,s2l1);
				qual2 = qual2.substr(0,s2l1);
			}

			//determine pairtype
			string pairtype = "ambi";
			if ( link1==1 && link2==1 ){
				pairtype = "expect";
				ebothcount++;
			} else
				pairtype = "ambi";
		
			if ( keepempty==2) { 
				if ((link1 == 0 && link2 == 1) ||
                	(link1 == 1 && link2 == 0) ) {
					pairtype = "expect";
					emptycount1++;
				} else if ( link1==0 && link2==0){
					pairtype = "expect";
					emptycount2++;
				}
			} else if ( keepempty==1){
				if ((link1 == 0 && link2 == 1) ||
                	(link1 == 1 && link2 == 0) ) {
					pairtype = "expect";
					emptycount1++;
				} else if ( link1==0 && link2==0){
					emptycount2++;
				}
			}
		
			if (pairtype == "expect")
				expectcount++;
        	else if (pairtype == "ambi")
				ambicount++;

			#pragma omp critical
			if ( (pairtype=="expect") &&  ( seq1.length() >= minlength ) && (seq2.length() >= minlength ))
			{
				validcount++;
				same1 << name1<<endl << seq1<<endl << "+" <<endl<< qual1<<endl;
				same2 << name2<<endl << seq2<<endl << "+" <<endl<< qual2<<endl;
			}

		}
		if (n_seqs1>0) cout << "Processed "<< n_seqs1<< " pair reads"<<endl<<flush;
		seqs1.clear();
		seqs2.clear();
	}
	gzclose(fp1);
	gzclose(fp2);
	same1.close();
	same2.close();
    cout << "Total PETs:\t" << linecount <<endl;
    cout << "Expect PETs:\t" << expectcount <<endl;
    cout << "Expect both PETs:\t" << ebothcount <<endl;
    cout << "Chim PETs:\t" << chimcount <<endl;
    //cout << "Ambi PETs: " << ambicount <<endl;
    cout << "1Empty PETs:\t" << emptycount1 <<endl;
    cout << "2Empty PETs:\t" << emptycount2 <<endl;
    cout << "Valid PETs:\t" << validcount <<endl<<flush;
    trimstat << "Total PETs:\t" << linecount <<endl;
    trimstat << "Expect PETs:\t" << expectcount <<endl;
    trimstat << "Expect both PETs:\t" << ebothcount <<endl;
    trimstat << "Chim PETs:\t" << chimcount <<endl;
    trimstat << "1Empty PETs:\t" << emptycount1 <<endl;
    trimstat << "2Empty PETs:\t" << emptycount2 <<endl;
    trimstat << "Valid PETs:\t" << validcount <<endl <<flush;
	trimstat.close();
	return 0;
}



vector< int > parseFastqgz(string fastqgz1, string fastqgz2,string basename,
              size_t minlength = 15,size_t maxlength = 300,
              int keepempty = 0, bool verbose = true,
              string linker1 = "GTTGGATAAG" , string linker2 = "GTTGGAATGT", int msk=0)
{

    // arguments
    gzFile file1 = gzopen(fastqgz1.c_str(),"r");
    gzFile file2 = gzopen(fastqgz2.c_str(),"r");
    ofstream same1 ( (basename + "_1.valid.fastq").c_str() );
    ofstream same2 ( (basename + "_2.valid.fastq").c_str() );
    //ofstream chim1 ( (basename + "_1.chim.fastq").c_str() );
    //ofstream chim2 ( (basename + "_2.chim.fastq").c_str() );
    
    // keep track of PET types
    int samecount = 0;
    int validcount = 0;
    int chimcount = 0;
    int ambicount = 0;
    int emptycount = 0;
    
    
    
    // define variables
    //string fqline1;
    //string fqline2;
    long long linecount = 0;
    int i = 0;
    vector< string > lines1;
    vector< string > lines2;
    int LENS = 500;
    string fqline1;
    string fqline2;
    char *buffer1 = new char[LENS];
    char *buffer2 = new char[LENS];
    //while (getline(file1, fqline1))
    while(0!=gzgets(file1,buffer1,LENS))
    {   
        fqline1 = buffer1;
        // read lines and increment counters
        //getline(file2, fqline2);
        gzgets(file2,buffer2,LENS);
        i++;
        linecount++;
        fqline2 = buffer2;
        // add lines to list
        lines1.push_back(fqline1);
        lines2.push_back(fqline2);
        // if list length is 4 perform operations, print, and clear lists
        if ( i == 4 )
        {
            // find the position of the linkers
            //int r1l1found = findlinker( lines1[1].c_str(), linker1.c_str(), msk);
            //int r1l2found = findlinker(lines1[1].c_str(), linker2.c_str(), msk);
            //int r2l1found = findlinker(lines2[1].c_str(), linker1.c_str(), msk);
            //int r2l2found = findlinker(lines2[1].c_str(), linker2.c_str(), msk);
            //cout << linecount/4 << ": "<<r1l1found << " " << r1l2found <<" " << r2l1found << " " << r2l2found <<endl;
            size_t r1l1found, r1l2found, r2l1found, r2l2found;
	    if (msk<1)
	    {
	    	r1l1found = lines1[1].find(linker1);
            	r1l2found = lines1[1].find(linker2);
            	r2l1found = lines2[1].find(linker1);
            	r2l2found = lines2[1].find(linker2);
            } else{
	        r1l1found = bitap_search(lines1[1], linker1, msk);
                r1l2found = bitap_search(lines1[1], linker2, msk);
                r2l1found = bitap_search(lines2[1], linker1, msk);
                r2l2found = bitap_search(lines2[1], linker2, msk);
            }
            
  	    // determine the linker type (0 = none, 3 = both)
            // read 1
            int r1linker = 0;
            if ( (r1l1found == npos) & (r1l2found == npos))
            {
                r1linker = 0;
            }
            else if ((r1l1found != npos) & (r1l2found == npos))
            {
                r1linker = 1;
                lines1[1] = lines1[1].substr(0,r1l1found).append("\n");
                lines1[3] = lines1[3].substr(0,r1l1found).append("\n");
            }
            else if ((r1l1found == npos) & (r1l2found != npos))
            {
                r1linker = 2;
                lines1[1] =  lines1[1].substr (0,r1l2found).append("\n");
                lines1[3] =  lines1[3].substr (0,r1l2found).append("\n");
            }
            else
            {
                r1linker = 3;
            }
            
            // read 2
            int r2linker = 0;
            if ((r2l1found == npos) & (r2l2found == npos))
            {
                r2linker = 0;
            }
            else if ((r2l1found != npos) & (r2l2found == npos))
            {
                r2linker = 1;
                lines2[1] =  lines2[1].substr (0,r2l1found).append("\n");
                lines2[3] =  lines2[3].substr (0,r2l1found).append("\n");
            }
            else if ((r2l1found == npos) & (r2l2found != npos))
            {
                r2linker = 2;
                lines2[1] =  lines2[1].substr (0,r2l2found).append("\n");
                lines2[3] =  lines2[3].substr (0,r2l2found).append("\n");
            }
            else
            {
                r2linker = 3;
            }
            
            // determine pairtype
            string pairtype = "ambi";
            if ((r1linker ==1 && r2linker==1 ) || (r1linker==2 && r2linker==2) )
            {
		pairtype = "same";
            }
            else if ( (r1linker==1 && r2linker==2) || (r1linker==2 && r2linker==1 ) )
            {
                pairtype = "chim";
            }
            else
            {
                pairtype = "ambi";   
            }

            if (keepempty == 2)
            {
                if ((r1linker == 0 && r2linker == 1) ||
                    (r1linker == 0 && r2linker == 2) ||
                    (r1linker == 1 && r2linker == 0) ||
                    (r1linker == 2 && r2linker == 0) ||
                    (r1linker == 0 && r2linker == 0))
                {
                    pairtype = "same";
		    emptycount++;
                }
	    } else if (keepempty == 1)
            {
                if ((r1linker == 0 && r2linker == 1) ||
                    (r1linker == 0 && r2linker == 2) ||
                    (r1linker == 1 && r2linker == 0) ||
                    (r1linker == 2 && r2linker == 0) )
                {
                    pairtype = "same";
 		    emptycount++;
                }
	    }
            // add to counters
            if (pairtype == "same")
            {
              samecount++;
            }
            if (pairtype == "chim")
            {
              chimcount++;
            }
            if (pairtype == "ambi")
            {
              ambicount++;
            }
            
            // determine if they pass the size requirements and print to output
            if ( (pairtype=="same") &&  ( lines1[1].length() >= minlength ) && (lines2[1].length() >= minlength ))
            {
		    validcount++;
                    same1 << lines1[0] << lines1[1] << lines1[2] << lines1[3];
                    same2 << lines2[0] << lines2[1] << lines2[2] << lines2[3];
            }
            
            // reset lines
            i = 0;
            lines1.clear();
            lines2.clear();

        }
        
        // num % 2 computes the remainder when num is divided by 2
        if ( linecount % 20000000 == 0 )
        {
            cout << linecount/4 <<endl;
	    cout << flush;
            
        }
    }
    
    // close streams
    same1.close();
    same2.close();
     
    // report results
    vector< int > parsingresults;
    parsingresults.push_back(samecount);
    parsingresults.push_back(chimcount);
    parsingresults.push_back(ambicount);
    cout << "Total PETs: " << linecount/4 <<endl;
    cout << "Same PETs: " << samecount <<endl;
    cout << "Chim PETs: " << chimcount <<endl;
    cout << "Ambi PETs: " << ambicount <<endl;
    cout << "Empty PETs: " << emptycount <<endl;
    cout << "Valid PETs: " << validcount <<endl;
    return parsingresults;
}



vector< int > parseFastqBridgegz(string fastqgz1, string fastqgz2,string basename,
              size_t minlength = 15, size_t maxlength = 300,
              int keepempty = 0, bool verbose = true,
              string linker1 = "ACGCGATATCTTATCTGACT" , string linker2 = "AGTCAGATAAGATATCGCGT", int msk=0)
{

    // arguments
    gzFile file1 = gzopen(fastqgz1.c_str(),"r");
    gzFile file2 = gzopen(fastqgz2.c_str(),"r");
    ofstream same1 ( (basename + "_1.valid.fastq").c_str() );
    ofstream same2 ( (basename + "_2.valid.fastq").c_str() );
    //ofstream chim1 ( (basename + "_1.chim.fastq").c_str() );
    //ofstream chim2 ( (basename + "_2.chim.fastq").c_str() );
    
    // keep track of PET types
    int samecount = 0;
    int validcount = 0;
    int chimcount = 0;
    int ambicount = 0;
    int emptycount = 0;
    
    // define variables
    //string fqline1;
    //string fqline2;
    long long linecount = 0;
    int i = 0;
    vector< string > lines1;
    vector< string > lines2;
    int LENS = 500;
    string fqline1;
    string fqline2;
    char *buffer1 = new char[LENS];
    char *buffer2 = new char[LENS];
    //while (getline(file1, fqline1))
    while(0!=gzgets(file1,buffer1,LENS))
    {   
        fqline1 = buffer1;
        // read lines and increment counters
        //getline(file2, fqline2);
        gzgets(file2,buffer2,LENS);
        i++;
        linecount++;
        fqline2 = buffer2;
        // add lines to list
        lines1.push_back(fqline1);
        lines2.push_back(fqline2);
        // if list length is 4 perform operations, print, and clear lists
        if ( i == 4 )
        {
            
            // find the position of the linkers
            //int r1l1found = findlinker( lines1[1].c_str(), linker1.c_str(), msk);
            //int r1l2found = findlinker(lines1[1].c_str(), linker2.c_str(), msk);
            //int r2l1found = findlinker(lines2[1].c_str(), linker1.c_str(), msk);
            //int r2l2found = findlinker(lines2[1].c_str(), linker2.c_str(), msk);
            //cout << linecount/4 << ": "<<r1l1found << " " << r1l2found <<" " << r2l1found << " " << r2l2found <<endl;
            size_t r1l1found, r1l2found, r2l1found, r2l2found;
	    if (msk <1)
	    {
	    	r1l1found = lines1[1].find(linker1);
            	r1l2found = lines1[1].find(linker2);
            	r2l1found = lines2[1].find(linker1);
            	r2l2found = lines2[1].find(linker2);
            } else{
	        r1l1found = bitap_search(lines1[1], linker1, msk);
                r1l2found = bitap_search(lines1[1], linker2, msk);
                r2l1found = bitap_search(lines2[1], linker1, msk);
                r2l2found = bitap_search(lines2[1], linker2, msk);
            }
  	    // determine the linker type (0 = none, 3 = both)
            // read 1
            int r1linker = 0;
            if ( (r1l1found == npos) & (r1l2found == npos))
            {
                r1linker = 0;
            }
            else if ((r1l1found != npos) & (r1l2found == npos))
            {
                r1linker = 1;
                lines1[1] =  lines1[1].substr (0,r1l1found).append("\n");
                lines1[3] =  lines1[3].substr (0,r1l1found).append("\n");
            }
            else if ((r1l1found == npos) & (r1l2found != npos))
            {
                r1linker = 2;
                lines1[1] =  lines1[1].substr (0,r1l2found).append("\n");
                lines1[3] =  lines1[3].substr (0,r1l2found).append("\n");
            }
            else
            {
                r1linker = 3;
            }
            
            // read 2
            int r2linker = 0;
            if ((r2l1found == npos) & (r2l2found == npos))
            {
                r2linker = 0;
            }
            else if ((r2l1found != npos) & (r2l2found == npos))
            {
                r2linker = 1;
                lines2[1] =  lines2[1].substr (0,r2l1found).append("\n");
                lines2[3] =  lines2[3].substr (0,r2l1found).append("\n");
            }
            else if ((r2l1found == npos) & (r2l2found != npos))
            {
                r2linker = 2;
                lines2[1] =  lines2[1].substr (0,r2l2found).append("\n");
                lines2[3] =  lines2[3].substr (0,r2l2found).append("\n");
            }
            else
            {
                r2linker = 3;
            }
            
            // determine pairtype
            string pairtype = "ambi";
            if ((r1linker ==1 && r2linker==2 ) || (r1linker==2 && r2linker==1) )
            {
		pairtype = "same";
            }
            else if ( (r1linker==1 && r2linker==1) || (r1linker==2 && r2linker==2 ) )
            {
                pairtype = "chim";
            }
            else
            {
                pairtype = "ambi";   
            }

            if (keepempty == 2)
            {
                if ((r1linker == 0 && r2linker == 1) ||
                    (r1linker == 0 && r2linker == 2) ||
                    (r1linker == 1 && r2linker == 0) ||
                    (r1linker == 2 && r2linker == 0) ||
                    (r1linker == 0 && r2linker == 0))
                {
                    pairtype = "same";
		    emptycount++;
                }
	    } else if (keepempty == 1)
            {
                if ((r1linker == 0 && r2linker == 1) ||
                    (r1linker == 0 && r2linker == 2) ||
                    (r1linker == 1 && r2linker == 0) ||
                    (r1linker == 2 && r2linker == 0) )
                {
                    pairtype = "same";
 		    emptycount++;
                }
	    }
            // add to counters
            if (pairtype == "same")
            {
              samecount++;
            }
            if (pairtype == "chim")
            {
              chimcount++;
            }
            if (pairtype == "ambi")
            {
              ambicount++;
            }
            
            // determine if they pass the size requirements and print to output
            if ( (pairtype=="same") && ( lines1[1].length() >= minlength ) && (lines2[1].length() >= minlength ))
            {
		validcount++;
                same1 << lines1[0] << lines1[1] << lines1[2] << lines1[3];
                same2 << lines2[0] << lines2[1] << lines2[2] << lines2[3];
            }
            // reset lines
            i = 0;
            lines1.clear();
            lines2.clear();
        }
        
        // num % 2 computes the remainder when num is divided by 2
        if ( linecount % 20000000 == 0 )
        {
            cout << linecount/4 <<endl;
	    cout << flush;
        }
    }
    
    // close streams
    same1.close();
    same2.close();
     
    // report results
    vector< int > parsingresults;
    parsingresults.push_back(samecount);
    parsingresults.push_back(chimcount);
    parsingresults.push_back(ambicount);
    cout << "Total PETs: " << linecount/4 <<endl;
    cout << "Same PETs: " << samecount <<endl;
    cout << "Chim PETs: " << chimcount <<endl;
    cout << "Ambi PETs: " << ambicount <<endl;
    cout << "Empty PETs: " << emptycount <<endl;
    cout << "Valid PETs: " << validcount <<endl;
    return parsingresults;
}


vector< int > parseFastqEnzymegz(string fastqgz1, string fastqgz2,string basename,
              size_t minlength = 15,size_t maxlength = 300,
              int keepempty = 0, bool verbose = true,
              string linker1 = "GATCGATC" , string linker2 = "GATCGATC", int msk=0)
{
    vector< int > parsingresults;
    if (linker1!=linker2){
	cout<< "Make sure linkerA is equal to linkerB in mode 2."<<endl;
	return parsingresults;
    }
    // arguments
    gzFile file1 = gzopen(fastqgz1.c_str(),"r");
    gzFile file2 = gzopen(fastqgz2.c_str(),"r");
    ofstream same1 ( (basename + "_1.valid.fastq").c_str() );
    ofstream same2 ( (basename + "_2.valid.fastq").c_str() );
    //ofstream chim1 ( (basename + "_1.chim.fastq").c_str() );
    //ofstream chim2 ( (basename + "_2.chim.fastq").c_str() );
    
    // keep track of PET types
    int samecount = 0;
    int validcount = 0;
    int chimcount = 0;
    int ambicount = 0;
    int emptycount = 0;
    
    // define variables
    //string fqline1;
    //string fqline2;
    long long linecount = 0;
    int i = 0;
    vector< string > lines1;
    vector< string > lines2;
    int LENS = 500;
    string fqline1;
    string fqline2;
    char *buffer1 = new char[LENS];
    char *buffer2 = new char[LENS];
    //while (getline(file1, fqline1))
    while(0!=gzgets(file1,buffer1,LENS))
    {   
        fqline1 = buffer1;
        // read lines and increment counters
        //getline(file2, fqline2);
        gzgets(file2,buffer2,LENS);
        i++;
        linecount++;
        fqline2 = buffer2;
        // add lines to list
        lines1.push_back(fqline1);
        lines2.push_back(fqline2);
        // if list length is 4 perform operations, print, and clear lists
        if ( i == 4 )
        {
            
            // find the position of the linkers
            //int r1l1found = findlinker( lines1[1].c_str(), linker1.c_str(), msk);
            //int r1l2found = findlinker(lines1[1].c_str(), linker2.c_str(), msk);
            //int r2l1found = findlinker(lines2[1].c_str(), linker1.c_str(), msk);
            //int r2l2found = findlinker(lines2[1].c_str(), linker2.c_str(), msk);
            //cout << linecount/4 << ": "<<r1l1found << " " << r1l2found <<" " << r2l1found << " " << r2l2found <<endl;
            size_t r1l1found, r2l1found;
	    if (msk<1)
	    {
	    	r1l1found = lines1[1].find(linker1);
            	r2l1found = lines2[1].find(linker1);
            } else{
	        r1l1found = bitap_search(lines1[1], linker1, msk);
                r2l1found = bitap_search(lines2[1], linker1, msk);
            }
            
  	    // determine the linker type (0 = none, 3 = both)
            // read 1
            int r1linker = 0;
            if ( r1l1found == npos )
            {
                r1linker = 0;
            }
            else
            {
                r1linker = 1;
                lines1[1] =  lines1[1].substr (0,r1l1found).append("\n");
                lines1[3] =  lines1[3].substr (0,r1l1found).append("\n");
            }
            
            // read 2
            int r2linker = 0;
            if ( r2l1found == npos )
            {
                r2linker = 0;
            }
            else
            {
                r2linker = 1;
                lines2[1] =  lines2[1].substr (0,r2l1found).append("\n");
                lines2[3] =  lines2[3].substr (0,r2l1found).append("\n");
            }
            
            // determine pairtype
            string pairtype = "ambi";
            if ( r1linker ==1 && r2linker==1 )
            {
		pairtype = "same";
            }
            else
            {
                pairtype = "ambi";   
            }

            if (keepempty == 2)
            {
                pairtype = "same";
		emptycount++;
	    } else if (keepempty == 1)
            {
                if ((r1linker == 0 && r2linker == 1) ||
                    (r1linker == 1 && r2linker == 0) )
                {
                    pairtype = "same";
 		    emptycount++;
                }
	    }
            // add to counters
            if (pairtype == "same")
            {
              samecount++;
            }
            if (pairtype == "ambi")
            {
              ambicount++;
            }
            
            // determine if they pass the size requirements and print to output
            if ( (pairtype=="same") && ( lines1[1].length() >= minlength ) && (lines2[1].length() >= minlength ))
            {
		validcount++;
                same1 << lines1[0] << lines1[1] << lines1[2] << lines1[3];
                same2 << lines2[0] << lines2[1] << lines2[2] << lines2[3];
            }
            
            // reset lines
            i = 0;
            lines1.clear();
            lines2.clear();
        }
        
        // num % 2 computes the remainder when num is divided by 2
        if ( linecount % 20000000 == 0 )
        {
            cout << linecount/4 <<endl;
	    cout << flush;
        }
    }
    
    // close streams
    same1.close();
    same2.close();
     
    // report results
    parsingresults.push_back(samecount);
    parsingresults.push_back(chimcount);
    parsingresults.push_back(ambicount);
    cout << "Total PETs: " << linecount/4 <<endl;
    cout << "Same PETs: " << samecount <<endl;
    cout << "Chim PETs: " << chimcount <<endl;
    cout << "Ambi PETs: " << ambicount <<endl;
    cout << "Empty PETs: " << emptycount <<endl;
    cout << "Valid PETs: " << validcount <<endl;
    return parsingresults;
}


typedef struct {
	int qual,start,stop;
	string name, seq, chrom,strand;
} align_t;

vector <align_t> read_sam(ifstream &file, int chunk_size, int *n_)
{
	int n=0;
	string line;
	vector <align_t> aligns;
	while (getline(file, line)) {
		align_t s;
        vector<string> es = string_split(line,"\t");
        s.name = es[0];
        int bitflag1 = StringToInt(es[1]);
        s.strand = get_strand(bitflag1);
        s.seq = es[9];
        s.chrom = es[2];
        s.start = StringToInt(es[3]) -1;
        s.stop = s.start + s.seq.length();
        s.qual = StringToInt(es[4]);
		aligns.push_back(s);
		n++;
		if (n >= chunk_size) break;
	}
	*n_ = n;
	return aligns;
}

void buildBedpe_p(string sam1, string sam2,string bedpefile, int mapq=30, int thread=1, int keepseq=0)
{
    omp_set_num_threads(thread);
    // arguments
    ifstream file1(sam1.c_str());
    ifstream file2(sam2.c_str());
    ofstream bedpefilestream ( bedpefile.c_str() );
	string statfile = bedpefile+".stat";
	ofstream statstream ( statfile.c_str() );
    // define variables
    string line1;
    string line2;
    long long linecount = 0;
	int lowmapN1 = 0;
	int unmapN1 = 0;
	int lowmapN2 = 0;
	int unmapN2 = 0;
	int outputN = 0;
    
    while (getline(file1, line1) )
    {
        // read lines and increment counter
        getline(file2, line2);
        if (line1.substr(0,1)!="@") 
			break;
	}
	int nalign1 = 1;
	int nalign2 = 1;

	while (nalign1>0){
		vector <align_t> aligns1 = read_sam(file1,2000000,&nalign1);
		vector <align_t> aligns2 = read_sam(file2,2000000,&nalign2);
		#pragma omp parallel for reduction(+:linecount) reduction(+:lowmapN1) reduction(+:unmapN1) reduction(+:lowmapN2) reduction(+:unmapN2) 
		for (int kk=0; kk<nalign1; kk++)
		{
			linecount++;
        	// split lines
			align_t rec1 = aligns1[kk];         
			align_t rec2 = aligns2[kk];         
        
			if (rec1.chrom == "*")
				unmapN1++;
			else if (rec1.qual<mapq)
				lowmapN1++;
			
			if (rec2.chrom == "*")
				unmapN2++;
			else if (rec2.qual<mapq)
				lowmapN2++;
        	// skip double stars
        	if ((rec1.chrom == "*") && (rec2.chrom == "*")){
          		continue;
			}


			// only use 3 kinds:
			// 60 60; * 60; 60 *;
			if ( !(rec1.qual>mapq && rec2.qual>mapq) && !(rec1.chrom=="*" && rec2.qual>mapq) && !(rec1.qual>mapq && rec2.chrom=="*") )
     		    continue;
       
			// check that read names match
        	if (rec1.name != rec2.name)
       		{
            	cout << "Error: read names of PET ends do not match";
            	exit (EXIT_FAILURE);
        	}
        
        	// determine which read goes first
        	bool reorder = false;
        	if ((rec1.chrom == rec1.chrom) && (rec1.start > rec2.start) )
            	reorder = true;
        	if ((rec1.chrom != rec2.chrom) && (rec1.chrom > rec2.chrom) )
            	reorder = true;
        	if ((rec1.chrom != rec2.chrom) && (rec1.chrom == "*") )
            	reorder = true;
        	if ((rec1.chrom != rec2.chrom) && (rec2.chrom == "*") )
            	reorder = false;
		
			// print out results
			string r1 = rec1.chrom+"\t"+IntToString(rec1.start)+"\t"+IntToString(rec1.stop)+"\t";
			string r2 = rec2.chrom+"\t"+IntToString(rec2.start)+"\t"+IntToString(rec2.stop)+"\t";
        	string line;
			if (reorder == false){
				if (keepseq==1)
            		line=r1+r2+rec1.name+"\t.\t"+rec1.strand+"\t"+rec2.strand+"\t"+rec1.seq+"\t"+rec2.seq+"\n";
				else
      				line=r1+r2+rec1.name+"\t.\t"+rec1.strand+"\t"+rec2.strand+"\n";
			}else{
				if (keepseq==1)
            		line=r2+r1+rec2.name+"\t.\t"+rec2.strand+"\t"+rec1.strand+"\t"+rec2.seq+"\t"+rec1.seq+"\n";
				else
            		line=r2+r1+rec2.name+"\t.\t"+rec2.strand+"\t"+rec1.strand+"\n";
			}
			#pragma omp critical
			{ bedpefilestream << line; outputN++;}
		}
		
		if (nalign1>0) cout << "Processed " << nalign1 << " PET" <<endl<<flush;
    }
    cout << "All\t" << linecount <<endl;
    cout << "Reads1 Low MAPQ\t" << lowmapN1 << "\t"<< setprecision(4) << (lowmapN1+0.0)/linecount*100<<"%"<<endl;
    cout << "Reads1 Unmapped\t" << unmapN1 << "\t"<< setprecision(4) << (unmapN1+0.0)/linecount*100<<"%"<<endl;
    cout << "Reads2 Low MAPQ\t" << lowmapN2 << "\t"<< setprecision(4) << (lowmapN2+0.0)/linecount*100<<"%"<<endl;
    cout << "Reads2 Unmapped\t" << unmapN2 << "\t"<< setprecision(4) << (unmapN2+0.0)/linecount*100<<"%"<<endl;
    cout << "Output PETs\t" << outputN << "\t"<< setprecision(4) << (outputN+0.0)/linecount*100<<"%"<<endl;
    statstream << "All\t" << linecount <<endl;
    statstream << "Reads1 Low MAPQ\t" << lowmapN1 <<endl;
    statstream << "Reads1 Unmapped\t" << unmapN1 <<endl;
    statstream << "Reads2 Low MAPQ\t" << lowmapN2 <<endl;
    statstream << "Reads2 Unmapped\t" << unmapN2 <<endl;
    statstream << "Output PETs\t" << outputN <<endl <<flush;
	statstream.close();
}




void buildBedpe(string sam1, string sam2,string bedpefile, int mapq)
{
    
    // arguments
    ifstream file1(sam1.c_str());
    ifstream file2(sam2.c_str());
    ofstream bedpefilestream ( bedpefile.c_str() );

    // define variables
    string line1;
    string line2;
    int linecount = 0;
    
    while (getline(file1, line1))
    {
        // read lines and increment counter
        getline(file2, line2);
        if (line1.substr(0,1)!="@")
	{
	
	linecount++;
        
        // split lines
        vector<string> e1 = string_split(line1,"\t");
        vector<string> e2 = string_split(line2,"\t");
        
        // get info for file 1
        string name1 = e1[0];
        //name1 = string_split(name1,"_")[0];
        name1 = string_split(name1," ")[0];
        //name1 = string_split(name1,"#")[0];
        int bitflag1 = StringToInt(e1[1]);
        string strand1 = get_strand(bitflag1);
        string sequence1 = e1[9];
        string chrom1 = e1[2];
        int start1 = StringToInt(e1[3]) -1;
        int qual1 = StringToInt(e1[4]);
        int stop1  = start1 + sequence1.length();
        
        // get info for file 2
        string name2 = e2[0];
        //name2 = string_split(name2,"_")[0];
        name2 = string_split(name2," ")[0];
        //name2 = string_split(name2,"#")[0];
        int bitflag2 = StringToInt(e2[1]);
        string strand2 = get_strand(bitflag2);
        string sequence2 = e2[9];
        string chrom2 = e2[2];
        int start2 = StringToInt(e2[3]) -1;
        int qual2 = StringToInt(e2[4]);
        int stop2  = start2 + sequence2.length();
        
        // skip double stars
        if ((chrom1 == "*") & (chrom2 == "*"))
        {
          continue;
        }
        
        // skip low qual
        if ((qual1<mapq) | (qual2<mapq))
        {
          continue;
        }
        
	// check that read names match
        if (name1 != name2)
        {
            cout << "Error: read names of PET ends do not match";
            break;
        }
        
        // determine which read goes first
        bool reorder = false;
        if ((chrom1 == chrom2) & (start1 > start2) )
        {
            reorder = true;
        }else if ((chrom1 != chrom2) & (chrom1 > chrom2) )
        {
            reorder = true;
        }else if ((chrom1 != chrom2) & (chrom1 == "*") )
        {
            reorder = true;
        }
	// print out results
        if (reorder == false)
        {
            vector<string> outputvector;
            outputvector.push_back(chrom1);
            outputvector.push_back(IntToString(start1));
            outputvector.push_back(IntToString(stop1));
            outputvector.push_back(chrom2);
            outputvector.push_back(IntToString(start2));
            outputvector.push_back(IntToString(stop2));
            outputvector.push_back(name1);
            outputvector.push_back(".");
            outputvector.push_back(strand1);
            outputvector.push_back(strand2);
            string outputstring = vector_join(outputvector,"\t");
            bedpefilestream << outputstring;
            bedpefilestream << "\n";
        }
        
        if (reorder == true)
        {
            vector<string> outputvector;
            outputvector.push_back(chrom2);
            outputvector.push_back(IntToString(start2));
            outputvector.push_back(IntToString(stop2));
            outputvector.push_back(chrom1);
            outputvector.push_back(IntToString(start1));
            outputvector.push_back(IntToString(stop1));
            outputvector.push_back(name1);
            outputvector.push_back(".");
            outputvector.push_back(strand2);
            outputvector.push_back(strand1);
            string outputstring = vector_join(outputvector,"\t");
            bedpefilestream << outputstring;
            bedpefilestream << "\n";
        }
	}
    }
}



void removeDups(string bedpein,string outnamebase,double distancesplit,int thread=1, int keepseq=0)
{
	// (1) split PETs by chrom and position 
    omp_set_num_threads(thread);
  	// keep track of output files
  	vector<string> outputvectorPETs;
  
  	// streams
 	ifstream file1(bedpein.c_str());
	string statfile = bedpein+".stat";
	ofstream statstream ( statfile.c_str(), ofstream::out | ofstream::app );
 	
	map<string, ofstream*> petsoutput;
  	string line;
  	while (getline(file1, line))
  	{
    	// split lines and determine bin
    	vector<string> currEall = string_split(line,"\t");
    	string chrom = currEall[0];
    	double pos = StringToInt(currEall[1]);
    	int bin = pos / distancesplit;
    	string binstring = NumberToString(bin);
    
    	// set output file name
    	string outname = outnamebase + "." + chrom + "." + binstring + ".bedpe";
    	// check if output file name exists (and make it fi neccesary)
    	if ( petsoutput.find(outname) == petsoutput.end() ) {  
      		petsoutput[outname] = new ofstream(outname.c_str());
      		outputvectorPETs.push_back( outname);  
    	}
    	// print to output file
    	*petsoutput[outname] << line <<endl;
	}
  	// close input stream
  	file1.close();
  
  	// close bedpe files streams
  	for (map<string, ofstream*>::iterator i = petsoutput.begin() ; i != petsoutput.end() ; i ++ ) {
    	i->second->close();
  	}    
  
  	// initialize counters
  	int nondups  = 0;
  	int duplines   = 0;
  	int interchrom = 0;
  	int intrachrom = 0;
  	int starchrom = 0;
  
  	// (2) read through each file and only print out non duplicates
  	// open input stream
  	string outputname = outnamebase;
  	ofstream finaloutput(outputname.c_str());
	int outputN = outputvectorPETs.size();
	#pragma omp parallel for reduction(+:nondups) reduction(+:duplines) reduction(+:interchrom) reduction(+:intrachrom) reduction(+:starchrom)
  	for (int i=0; i<outputN; i ++ ) {
    	// make new map
    	map<string, int> PETmap;
    
    	// open input stream
    	ifstream file1(outputvectorPETs[i].c_str());
    	//ifstream file1(i->c_str());
    
    	string line;
		//#pragma omp critical
   		while (getline(file1, line))
   		{
      		// split lines and determine bin
      		vector<string> currEall = string_split(line,"\t");
      		string chrom1  = currEall[0];
      		string pos1  = currEall[1];
      		string pos11  = currEall[2];
      		string chrom2 = currEall[3];
      		string pos2 = currEall[4];
      		string pos22 = currEall[5];
      		string uniqcode = chrom1+"_"+pos1+"_"+pos11+"_"+chrom2+"_"+pos2+"_"+pos22;

      		// check if read has been seen before
      		if ( PETmap.find(uniqcode) == PETmap.end() ) {  
        		PETmap[uniqcode] = 0;
      		}
      		PETmap[uniqcode]++;
      
      		if (PETmap[uniqcode] > 1)
      		{
        		duplines++;
      		}
      
      		// if it is the first instance print it to the output file
      		if (PETmap[uniqcode] == 1)
      		{
        		nondups++;
        		if (chrom1 == chrom2)
          			intrachrom++;
				else if (chrom1=="*" || chrom2=="*")
					starchrom++;
				else
          			interchrom++;
				line.append("\n");
		 		#pragma omp critical
        		finaloutput << line;
        	}
      	}
    	// close input stream
    	file1.close();
    	//remove(i->c_str());
    	remove(outputvectorPETs[i].c_str());
  	}
  	// close input stream
  	finaloutput.close();
  
  	// report results
    cout << "Duplicates\t" << duplines <<endl;
    cout << "Uniques\t" << nondups <<endl;
    cout << "Intra\t" << intrachrom <<endl;
    cout << "Inter\t" << interchrom <<endl;
    cout << "OneEndMapped\t" << starchrom <<endl <<flush;
    statstream << "Duplicates\t" << duplines <<endl;
    statstream << "Uniques\t" << nondups <<endl;
    statstream << "Intra\t" << intrachrom <<endl;
    statstream << "Inter\t" << interchrom <<endl;
    statstream << "OneEndMapped\t" << starchrom <<endl <<flush;
	statstream.close();
}


void buildTagAlign(string bedpefile, string TagAlignfile) {
    
    // establish streams
    ifstream infile (bedpefile.c_str());
    ofstream outfile ( TagAlignfile.c_str() );

   int i = 0;
   string line;
   while (getline(infile, line))
    {
        // increment counters
        i++;
        // split line by tab
         vector<string> e = string_split(line,"\t");
        
        // reverse strands
        string newstrand1 = "-";
        if (e[8] == "-")
        {
          newstrand1 = "+";
        }
        e[8] = newstrand1;
        
        string newstrand2 = "-";
        if (e[9] == "-")
        {
          newstrand2 = "+";
        }
        e[9] = newstrand2;
        
        // print to output
        if (e[0] != "*")
        {
          vector<string> outputvector;
          outputvector.push_back(e[0]);
          outputvector.push_back(e[1]);
          outputvector.push_back(e[2]);
          outputvector.push_back(e[6]);
          outputvector.push_back(".");
          outputvector.push_back(e[8]);
          string outputstring = vector_join(outputvector,"\t");
          outfile << outputstring<<endl;
        }
        
        if (e[3] != "*")
        {
          vector<string> outputvector2;
          outputvector2.push_back(e[3]);
          outputvector2.push_back(e[4]);
          outputvector2.push_back(e[5]);
          outputvector2.push_back(e[6]);
          outputvector2.push_back(".");
          outputvector2.push_back(e[9]);
          string outputstring2 = vector_join(outputvector2,"\t");
          outfile << outputstring2 <<endl;
        }
    }
    // close streams
    infile.close();
    outfile.close();
}



void DeterminePeakDepths(string temppeakoverlap,string peaksfileslopdepth)
{
  // streams
  ifstream fileIN(temppeakoverlap.c_str());
  ofstream fileOUT(peaksfileslopdepth.c_str());
  
  // make map of peaks
  map<string, int > peaksmap;
  
  // read in file line by line
  string line;
  while (getline(fileIN, line))
  {
    // split lines and determine bin
    vector<string> currEall = string_split(line,"\t");
    string peakchrom = currEall[6];
    string peakstart = currEall[7];
    string peakend   = currEall[8];
    string peakname  = currEall[9];
  
    // peak name
    string longname = peakchrom + "\t" + peakstart + "\t" + peakend + "\t" + peakname;
  
    // add peak to map
    if ( peaksmap.find(longname) == peaksmap.end() ) {  
      peaksmap[longname] = 0;
    }
    
    // increment count
    peaksmap[longname]++;
    
  }
  
  // close input file
  fileIN.close();
  
  // print out info
  for (map<string, int >::const_iterator longname = peaksmap.begin() ; longname != peaksmap.end() ; longname ++ ){
    fileOUT << longname->first + "\t" + NumberToString(longname->second) + "\t.";
    fileOUT << "\n";
  }
  
  // close output file
  fileOUT.close();

}

void bedpe2Interaction(string bedpefile, string peaksfile, string statfile)
{
  // streams
  ifstream fileIN(bedpefile.c_str());
  string fintra = peaksfile+".intra.bedpe";
  string finter = peaksfile+".inter.bedpe";
  ofstream fileOUT1(fintra.c_str());
  ofstream fileOUT2(finter.c_str());
  
  // make map of peaks
  map<string, int > peaksmap1;
  map<string, int > peaksmap2;
  
  // read in file line by line
  vector<string> lines;
  string line;
  int selfligationN = 0;
  //int singletonN = 0;
  //int interligationN = 0;
  int i= 0;
  while (getline(fileIN, line))
  {
    lines.push_back(line);
    i++;
    if (i==2)
    {
		i = 0;
	// split lines and determine bin
    	vector<string> currEalla = string_split(lines[0],"\t");
    	vector<string> currEallb = string_split(lines[1],"\t");
		string peakchrom1, peakstart1, peakend1, peakname1, peakdep1;
		string peakchrom2, peakstart2, peakend2, peakname2, peakdep2;
		if (currEalla.size()>16){ // for phased bedpe
    		peakchrom1 = currEalla[12];
    		peakstart1 = currEalla[13];
    		peakend1   = currEalla[14];
    		peakname1  = currEalla[15];
    		peakdep1  = currEalla[16];
    		
			peakchrom2 = currEallb[12];
    		peakstart2 = currEallb[13];
    		peakend2   = currEallb[14];
    		peakname2  = currEallb[15];
    		peakdep2  = currEallb[16];
		} else{
    		peakchrom1 = currEalla[10];
    		peakstart1 = currEalla[11];
    		peakend1   = currEalla[12];
    		peakname1  = currEalla[13];
    		peakdep1  = currEalla[14];
    		
			peakchrom2 = currEallb[10];
    		peakstart2 = currEallb[11];
    		peakend2   = currEallb[12];
    		peakname2  = currEallb[13];
    		peakdep2  = currEallb[14];
		}
		
		if (currEalla[1] != currEallb[1]){
			//cout << "Pair overlap ERROR!!!"<<endl;
			//cout << peakname1 << "---" << peakname2<<endl;
			string tline = lines[1];
			lines.clear();
    		lines.push_back(tline);
			i = 1;
			continue;
			//exit (EXIT_FAILURE);
		}

    
    	if (peakname1==peakname2) {
			selfligationN++;
			lines.clear();
			continue;
    	}
	
    	// interaction name
    	string longname = peakchrom1 + "\t" + peakstart1 + "\t" + peakend1 +"\t"+peakchrom2 + "\t" + peakstart2 + "\t" + peakend2 + "\t" + peakname1+"\t"+peakname2+ "\t"+ peakdep1 + "\t"+ peakdep2;
    	// add peak to map
    	if ( peakchrom1==peakchrom2 && peaksmap1.find(longname) == peaksmap1.end() ) {  
      		peaksmap1[longname] = 1;
   	 	} else if ( peakchrom1==peakchrom2)
    		peaksmap1[longname]++;
		else if ( peaksmap2.find(longname) == peaksmap2.end() ) {  
      		peaksmap2[longname] = 1;
   	 	} else {
			peaksmap2[longname]++;
		}// increment count
 		lines.clear();
     }
  }
  
  // close input file
  fileIN.close();
  int intraN=0;
  int interN=0;
  int intraloopN=0;
  int interloopN=0;

  // print out info
  for (map<string, int >::const_iterator longname = peaksmap1.begin() ; longname != peaksmap1.end() ; longname ++ ){
    fileOUT1 << longname->first + "\t" + NumberToString(longname->second) <<endl;
	intraN = intraN + longname->second;
	intraloopN++;
  }
  for (map<string, int >::const_iterator longname = peaksmap2.begin() ; longname != peaksmap2.end() ; longname ++ ){
    fileOUT2 << longname->first + "\t" + NumberToString(longname->second) << endl;
	interN = interN + longname->second;
	interloopN++;
  }
  
  // close output file
  fileOUT1.close();
  fileOUT2.close();
  ofstream statout( statfile.c_str(), ofstream::out | ofstream::app );
  statout << "PETs in the same peak\t" << selfligationN <<endl;
  statout << "Intra PETs bewteen two peaks\t" << intraN <<endl;
  statout << "Inter PETs bewteen two peaks\t" << interN <<endl;
  statout.close();
  cout << "PETs in the same peak: " << selfligationN <<endl;
  cout << "Intra PETs bewteen two peaks: " << intraN <<endl;
  cout << "Inter PETs bewteen two peaks: " << interN <<endl;
  cout << "Intra loops: " << intraloopN <<endl;
  cout << "Inter loops: " << interloopN <<endl<<endl<<flush;
}

void bedpe2Phased(string bedpefile, string output)
{
	// streams
  	ifstream fileIN(bedpefile.c_str());
  	ofstream fileOUT(output.c_str());
  
  	// make map of peaks
  	map<string, string > map1;
  
  	// read in file line by line
  	vector<string> lines;
  	string line;
	int N = 0, Nc = 0;
  	while (getline(fileIN, line))
  	{
		// split lines, expected 17=12+5 columns
    	vector<string> currEalla = string_split(line,"\t");
    	string chrom1 = currEalla[0];
    	int start1 = StringToInt(currEalla[1]);
    	int end1 = StringToInt(currEalla[2]);
    	string chrom2 = currEalla[3];
    	int start2 = StringToInt(currEalla[4]);
    	int end2 = StringToInt(currEalla[5]);
    	string seq1 = currEalla[10];
    	string seq2 = currEalla[11];
    	
		//skip "chrom *"
		if (chrom1=="*" || chrom2=="*") 
			continue;

		string vcfchrom = currEalla[12];
    	int vcfstart = StringToInt(currEalla[13]);
    	//int vcfend = StringToInt(currEalla[14]);
    	string vcfp   = currEalla[15];
    	string vcfm  = currEalla[16];
		
		int vflag1 = 0;
		int vflag2 = 0;
    	if (chrom1==vcfchrom && vcfstart>=start1 && vcfstart<end1) {
			int id = vcfstart-start1;
			string snp = seq1.substr(id,1);
			if (snp==vcfp)
				vflag1 = 1;
    		else if (snp==vcfm)
				vflag1 = 2;
			else
				vflag1 = 0;
		}
		if (chrom2==vcfchrom && vcfstart>=start2 && vcfstart<end2) {
			int id = vcfstart-start2;
			string snp = seq2.substr(id,1);
			if (snp==vcfp)
				vflag2 = 1;
    		else if (snp==vcfm)
				vflag2 = 2;
			else 
				vflag2 = 0;
		}
		
		if (vflag1+vflag2==3)
			continue;
			
    	// interaction name
    	//string longname = chrom1 + "\t" + IntToString(start1) + "\t" + IntToString(end1) +"\t"+chrom2 + "\t" + IntToString(start2) + "\t" + IntToString(end2) + "\t" + seq1+"\t"+seq2+ "\t"+ vcfchrom+ "\t"+ IntToString(vcfstart)+"\t"+ IntToString(vcfend) + "\t"+ vcfp + "\t"+ vcfm;
    	string longname = chrom1 + "\t" + IntToString(start1) + "\t" + IntToString(end1) +"\t"+chrom2 + "\t" + IntToString(start2) + "\t" + IntToString(end2);
    	// add peak to map
		if (vflag1+vflag2 > 0) {
			if ( map1.find(longname) == map1.end() ) {// not in map
      			map1[longname] = IntToString(vflag1)+"\t"+IntToString(vflag2);
			} else if (map1[longname] != "conflict"){
				vector <string> vals = string_split(map1[longname],"\t");
				int v1 = StringToInt(vals[0]);
				int v2 = StringToInt(vals[1]);
				if( v1+vflag1==3 || v2+vflag2==3 ) {
					map1[longname] = "conflict";
				} else {
					if (v1<vflag1) v1=vflag1;
					if (v2<vflag2) v2=vflag2;
					map1[longname] = IntToString(v1)+"\t"+IntToString(v2);
				}
			}
		}
	} //end while
  	// close input file
  	fileIN.close();

  	// print out info
  	for ( map <string,string >::const_iterator i = map1.begin() ; i != map1.end() ; i++ ){
    	if ( i->second != "conflict" ) {
			fileOUT << i->first + "\t" + i->second <<endl;
			N++;
		}else Nc++;
  	}
  
  	// close output file
  	fileOUT.close();
  	cout << "Phased PETs: " << N <<endl;
  	cout << "Conflict PETs: " << Nc <<endl;
}

