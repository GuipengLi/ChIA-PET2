#include <iostream>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "ChIA-PET2.h"
using namespace std;

/*funcion that show the help information*/
void usage(char *s)
{
  cout<<"Usage:   "<<s<<" sam1 sam2 out mapq_cutoff thread keepseq_flag"<<endl;
}


int main (int argc,char *argv[])
{
  // parameters
  string sam1;
  string sam2;
  string outbedpe;
  int mapq = 30;
  int thread = 1;
  int keepseq_flag = 0;
  /*if the program is ran witout options ,it will show the usgage and exit*/
  if(argc == 1)
  {
    usage(argv[0]);
    exit(1);
  }
  /*use function getopt to get the arguments with option."hu:p:s:v" indicate 
  that option h,v are the options without arguments while u,p,s are the
  options with arguments*/
  if (argc == 4) {
	sam1 = argv[1];
	sam2 = argv[2];
	outbedpe = argv[3];
  }else if (argc ==5){
	sam1 = argv[1];
	sam2 = argv[2];
	outbedpe = argv[3];
	mapq = atoi(argv[4]);
  }else if (argc ==6){
	sam1 = argv[1];
	sam2 = argv[2];
	outbedpe = argv[3];
	mapq = atoi(argv[4]);
	thread = atoi(argv[5]);
  }else if (argc ==7){
	sam1 = argv[1];
	sam2 = argv[2];
	outbedpe = argv[3];
	mapq = atoi(argv[4]);
	thread = atoi(argv[5]);
	keepseq_flag = atoi(argv[6]);
  }else{
     cout <<"Wrong parameters."<<endl;
     usage(argv[0]);
     exit(1);
  }

  //start
  cout <<"Running buildbedpe..."<<endl;
  buildBedpe_p(sam1, sam2, outbedpe, mapq,thread, keepseq_flag);
  return 0;
}


