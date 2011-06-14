/* A program to test cross-correlations. Written in OpenMP parallel C++. Scott Michael.*/
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stack>

const int minrows=30;
const int selectcell=23;

template<typename T, int size>
int GetArrLength(T(&)[size]){return size;}

using namespace std;
struct sample {
  float plot;
  float ba;
  int sp;
  int tree;
  int seed;
  int cell;
};

void sift_down(sample &, const int, const int);
void hpsort(sample &);


int main(int argc, char *argv[])
{

  sample *inputdata, *subdata;
  string filename,inputstring,substring;
  ifstream input;
  int c=0,b,i,inputsize,unique,subdatasize;
  stack<int> tempstack;
  int *spuni;

  filename = "sample.input";
  
  //START GETFILE FUNCTION
  FILE *f=fopen(filename.c_str(),"rb");
  while((b=fgetc(f))!=EOF) c+=(b==10)?1:0;
  fclose(f);
  
  inputsize = c;
  inputdata = new sample[inputsize];
  c=0;
    
  input.open(filename.c_str());
  getline(input, inputstring);
  while(getline(input, inputstring)){
    istringstream iss(inputstring);
    iss >> inputdata[c].plot;
    iss >> inputdata[c].sp;
    iss >> inputdata[c].tree;
    iss >> inputdata[c].ba;
    iss >> inputdata[c].seed;
    iss >> inputdata[c].cell;
    c++;
  }
  
  i=0;
  for(c=0;c<inputsize-1;c++)
    {
      if(inputdata[c].cell == selectcell) i++;
    }
  subdata = new sample[i];
  i=0;
  for(c=0;c<inputsize-1;c++){
    if(inputdata[c].cell == selectcell){
      subdata[i] = inputdata[c];
      i++;
    }
  }
  
  delete [] inputdata;

  subdatasize = i;
  hpsort(subdata);
  tempstack.push(subdata[0].sp);
  unique = 1;
  for(i=1;i<subdatasize;i++)
    {
      if(subdata[i].sp!=subdata[i-1].sp) unique++;
      tempstack.push(subdata[i].sp);
    }
  
  spuni = new int[unique];
  
  for(i=0;i<unique;i++)
    {
      spuni[i] = tempstack.top();
      tempstack.pop();
    }

  for(i=0;i<unique;i++)
    {
      

  
      

  

  //FILE GOTTEN

  


  

  return 0;
}
  
void sift_down(sample &tosort, const int l, const int r)
{
  int j, jold;
  sample a;
  
  a = tosort[l];
  jold=l;
  j=l+1;
  while (j <= r){
    if (j < r && tosort[j].sp < tosort[j+1].sp) j++;
    if (a.sp >= tosort[j].sp) break;
    tosort[jold]=tosort[j];
    jold=j;
    j=2*j+1;
  }
  tosort[jold]=a;
}

void hpsort(sample &tosort)
{
  int i;
  int n = GetArrLength(tosort);
  for (i=n/2-1; i>=0; i--) sift_down(tosort, i n-1);
  for (i=n-1; i>0, i--){
    SWAP(tosort[0],tosort[i]);
    sift_down(tosort,0,i-1);
  }
}
  
