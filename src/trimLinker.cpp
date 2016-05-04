#include <iostream>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "ChIA-PET2.h"
using namespace std;

/*funcion that show the help information*/
void usage(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-t  threads"<<endl;
  cout<<"         "<<"-e  number of mismatch allowed in linker, default=0"<<endl;
  cout<<"         "<<"-k  keep empty: 0, 1, 2 "<<endl;
  cout<<"         "<<"-l  min length of timmed reads: default=15"<<endl;
  cout<<"         "<<"-o  output directory"<<endl;
  cout<<"         "<<"-m  mode: 0 or 1, A/B linkers(0) or bridge linker(1)"<<endl;
  cout<<"         "<<"-A  linkerA"<<endl;
  cout<<"         "<<"-B  linkerB"<<endl;
  cout<<"         "<<"-n  output name prefix"<<endl;
  cout<<"example: "<<s<<" -t threads -o outdir -m mode -n name -k 0 -A linkerA -B likerB fq1 fq2"<<endl;
}


int main (int argc,char *argv[])
{
  char tmp;
  // parameters
  int THREAD = 1;
  int keepempty = 0;
  size_t errornum = 0;
  string linkerA ="GTTGGATAAG";
  string linkerB ="GTTGGAATGT";
  string OUTDIR = "output";
  string OUTNAME = "out";
  string fq1;
  string fq2;
  int mode = 0;
  int length = 15;
  /*if the program is ran witout options ,it will show the usgage and exit*/
  if(argc == 1)
  {
    usage(argv[0]);
    exit(1);
  }
  /*use function getopt to get the arguments with option."hu:p:s:v" indicate 
  that option h,v are the options without arguments while u,p,s are the
  options with arguments*/
  while((tmp=getopt(argc,argv,"hA:B:t:m:e:k:l:o:n:"))!=-1)
  {
    switch(tmp)
    {
      /*option h show the help infomation*/
      case 'h':
        usage(argv[0]);
        break;
      /*option u present the username*/
      case 't':
	THREAD = atoi(optarg);
        cout<<"thread is "<< THREAD<<endl;
        break;
      case 'm':
	mode = atoi(optarg);
        cout<<"mode is "<< mode<<endl;
        break;
      case 'k':
	keepempty = atoi(optarg);
        cout<<"keepempty is "<< keepempty<<endl;
        break;
      case 'l':
	length = atoi(optarg);
        cout<<"min length of trimmed read is "<< length<<endl;
        break;
      case 'e':
	errornum= atoi(optarg);
        cout<<"Error allowed is "<< errornum<<endl;
        break;
      case 'A':
	linkerA = optarg;
        cout<<"linkerA is "<< linkerA <<endl;
        break;
      case 'B':
	linkerB = optarg;
        cout<<"linkerB is "<< linkerB <<endl;
        break;
      case 'o':
	OUTDIR = optarg;
        cout<<"Output dir is "<< OUTDIR <<endl;
        break;
      case 'n':
	OUTNAME = optarg;
        cout<<"Output name is "<<OUTNAME <<endl;
	break;
      /*invail input will get the heil infomation*/
      default:
        usage(argv[0]);
      break;
    }
  }
  if (argc - optind ==2) {
	fq1 = argv[optind];
	fq2 = argv[optind+1];
	cout << "Reads 1: "<< fq1 <<endl;
	cout << "Reads 2: "<< fq2 <<endl;
  }else{
     cout <<"Wrong parameters."<<endl;
     usage(argv[0]);
     exit(1);
  }

  //start
  string basename = OUTDIR+"/"+OUTNAME;
  //if (mode==0)
  //    vector <int> trimresults = parseFastqgz(fq1, fq2, basename, length, 300, keepempty, true,linkerA, linkerB, errornum);
  //else if (mode==1)
  //    vector <int> trimresults = parseFastqBridgegz(fq1, fq2, basename, length, 300, keepempty, true,linkerA, linkerB, errornum);
  //else if (mode==2)
  //    vector <int> trimresults = parseFastqEnzymegz(fq1, fq2, basename, length, 300, keepempty, true,linkerA, linkerB, errornum);
  if (mode==0)
      parseFastqgz_p(THREAD, fq1, fq2, basename, length, keepempty, linkerA, linkerB, errornum);
  else if (mode==1)
      parseFastqBridgegz_p(THREAD,fq1, fq2, basename, length, keepempty,linkerA, linkerB, errornum);
  else if (mode==2)
      parseFastqEnzymegz_p(THREAD,fq1, fq2, basename, length,keepempty, linkerA, linkerB, errornum);
  else 
      {cout << "Choose mode: 0 or 1 or 2"<<endl;exit(1);}

  //cout << trimresults[0] <<endl;
  return 0;
}


