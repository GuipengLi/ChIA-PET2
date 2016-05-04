#include <iostream>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "ChIA-PET2.h"
using namespace std;

/*funcion that show the help information*/
void usage(char *s)
{
  cout<<"Usage:   "<<s<<" tagin depthout"<<endl;
}


int main (int argc,char *argv[])
{
  // parameters
  string tagin;
  string depthout;
  /*if the program is ran witout options ,it will show the usgage and exit*/
  if(argc == 1)
  {
    usage(argv[0]);
    exit(1);
  }
  /*use function getopt to get the arguments with option."hu:p:s:v" indicate 
  that option h,v are the options without arguments while u,p,s are the
  options with arguments*/
  if (argc == 3) {
	tagin = argv[1];
	depthout = argv[2];
  }else{
     cout <<"Not enough parameters."<<endl;
     usage(argv[0]);
     exit(1);
  }

  //start
  cout <<"Running tag2depth..."<<endl;
  DeterminePeakDepths(tagin, depthout);  
  return 0;
}


