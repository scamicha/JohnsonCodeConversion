/* A program to test cross-correlations. Written in C++ by Scott Michael. Converted from an R
 code provided by Daniel Johnson. The code takes as input a file containing the data plot, sp,
 tree, ba, seed, and cell in column format. Input and output files are specified in the
 source. In addition, the number of random iterations, the cell ofinterest and the minimum
 number of sp's to consider are specified in the source.

 compile this code by: g++ -O3 -o <executable name> sample.cpp

 The executable takes two required arguments and one optional argument and will only write
 errors to STDOUT. These errors can be thrown by the random number generator. The arguments
 are the input file, the output file, and optionally the number of iterations. The code is
 single precision and will only generate a few million random numbers. i.e. Don't set 
 iter > 10e6 without a rewrite of runif and ran1.*/

//header includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stack>
#include <iomanip>
#include <math.h>

using namespace std;

// These are parameters that can be changed by the user. Changing them requires the code to be
// recompiled. Iter can be entered at the command line if desired without recompilation. If it
// is not present then it defaults to the value here.
int iter = 5000;
const int minrows=30;
// End user modifiable parameters.

long int idum = 0;

// Struct to hold input data.
struct sample {
  float plot;
  float ba;
  int sp;
  int tree;
  int seed;
  int cell;
};

// Functions at the end of the source.
void sift_down(sample [], const int, const int);  // Used to sort input data by sp 
void hpsort(sample []);   // Function to sort data by sp
int *randgen(int);        // Generates an array of random integers
double ran1();            // Random number generator
float runif(float,float); // Random number generator


int main(int argc, char *argv[])
{

  sample *inputdata=NULL;
  sample *subdata=NULL;
  string inputstring,substring,filename,outfile;
  if(argc != 3 && argc != 4){
    cout << "Incorrect number of inputs. Aborting." << endl;
    return 1;
  }
  ifstream input;
  ofstream output;
  int c=0,b,i,inputsize,unique,subdatasize,unistart,uniend;
  int treesize,treesize1,x,y,j,k,l,yr,treesum,seedsum;
  stack<int> tempstack;
  int *spuni=NULL;
  int *seed=NULL;
  int *tree=NULL;
  int *randsel=NULL;
  int *seedran=NULL;
  int *spused=NULL;
  float q[4],eq[4],xdum[2],ydum[2];
  float *dbks=NULL;
  float *pvalue=NULL;
  float *dbksran=NULL;

  input.open(argv[1],ios::in);
  output.open(argv[2],ios::out);
  if (argc == 4) iter = atoi(argv[3]);

  c=0;
  while(!input.eof()){
    getline(input,inputstring);
    c++;
  }
  c--;
  input.clear();
  input.seekg(0, ios::beg);  
  
  inputsize = c;
  inputdata = new sample[inputsize];
  c=0;
 
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
      subdata[i] = inputdata[c];
      i++;
  }
  
  delete [] inputdata;
  inputdata=NULL;
  input.close();

  subdatasize = i;
  hpsort(subdata);
  tempstack.push(subdata[0].sp);
  unique = 1;
  for(i=1;i<subdatasize;i++)
    {
      if(subdata[i].sp!=subdata[i-1].sp){
	unique++;
	tempstack.push(subdata[i].sp);
      }
    }
  
  spuni   = new int[unique];
  dbks    = new float[unique];
  pvalue  = new float[unique];
  spused  = new int[unique];
  dbksran = new float[iter]; 
  
  for(i=0;i<unique;i++)
    {
      spuni[i] = tempstack.top();
      tempstack.pop();
      dbks[i] = 0.0;
      pvalue[i] = 0.0;
      spused[i] = 0;
    }
 
  for(i=0;i<unique/2;i++)
    swap(spuni[i],spuni[unique-1-i]);

  c=0;

  for(i=0;i<unique;i++){
      for(j=0;j<iter;j++)
	dbksran[j] = 0.0;
      unistart=c;
      while(c<subdatasize-1)
	if (subdata[c].sp == subdata[c+1].sp)
	  c++;
	else
	  break;
      c++;
      uniend=c;
      treesize = uniend-unistart;
      if (treesize < minrows)
	continue;
      seedran = new int[treesize];
      seed    = new int[treesize];
      tree    = new int[treesize];
      
      b=0;
      for(j=unistart;j<uniend;j++)
	{
	  seed[b] = subdata[j].seed;
	  tree[b] = subdata[j].tree;
	  b++;
	}
      treesum = 0;
      seedsum = 0;
      for(j=0;j<treesize;j++){
	treesum += tree[j];
	seedsum += seed[j];
      }
      
      if (treesum == 0 || seedsum == 0){
	delete [] seedran;
	seedran=NULL;
	delete [] seed;
	seed=NULL;
	delete [] tree;
	tree=NULL;
	continue;
      }
      spused[i] = 1;
      
      for(j=0;j<treesize;j++)
	{
	  x = tree[j];
	  y = seed[j];
	  for(k=0;k<4;k++){
	    q[k]=0.0;
	    eq[k]=0.0;
	  }
	  for(k=0;k<2;k++){
	    xdum[k]=0.0;
	    ydum[k]=0.0;
	  }

	  for(k=0;k<treesize;k++)
	    {
	      if(tree[k]>x && seed[k]>=y) q[0]+=1.0;
	      if(tree[k]<=x && seed[k]>=y) q[1]+=1.0;
	      if(tree[k]<=x && seed[k]<y) q[2]+=1.0;
	      if(tree[k]>x && seed[k]<y) q[3]+=1.0;
	      if(tree[k]>x){
		xdum[0]+=1.0;
	      }else xdum[1]+=1.0;
	      if(seed[k]>=y){
		ydum[0]+=1.0;
	      }else ydum[1]+=1.0;
	    }
	  eq[0] = xdum[0]*ydum[0]/(treesize*treesize);
	  eq[1] = xdum[1]*ydum[0]/(treesize*treesize);
	  eq[2] = xdum[1]*ydum[1]/(treesize*treesize);
	  eq[3] = xdum[0]*ydum[1]/(treesize*treesize);
	  for(k=0;k<4;k++){
	    q[k] /= treesize;
	    if(fabs(q[k]-eq[k]) > dbks[i]) 
	      dbks[i] = fabs(q[k]-eq[k]);
	  }
	}    
    
      for(l=0;l<iter;l++){
	randsel = randgen(treesize);
	for(j=0;j<treesize;j++)
	  seedran[j] = seed[randsel[j]];
	for(j=0;j<treesize;j++)
	  {
	    x = tree[j];
	    yr = seedran[j];
	    for(k=0;k<4;k++){
	      q[k]=0.0;
	      eq[k]=0.0;
	    }
	    for(k=0;k<2;k++){
	      xdum[k]=0.0;
	      ydum[k]=0.0;
	    }
	    
	    for(k=0;k<treesize;k++)
	      {
		if(tree[k]>x && seedran[k]>=yr) q[0]+=1.0;
		if(tree[k]<=x && seedran[k]>=yr) q[1]+=1.0;
		if(tree[k]<=x && seedran[k]<yr) q[2]+=1.0;
		if(tree[k]>x && seedran[k]<yr) q[3]+=1.0;
		if(tree[k]>x){
		  xdum[0]+=1.0;
		}else xdum[1]+=1.0;
		if(seedran[k]>=yr){
		  ydum[0]+=1.0;
		}else ydum[1]+=1.0;
	      }
	    eq[0] = xdum[0]*ydum[0]/(treesize*treesize);
	    eq[1] = xdum[1]*ydum[0]/(treesize*treesize);
	    eq[2] = xdum[1]*ydum[1]/(treesize*treesize);
	    eq[3] = xdum[0]*ydum[1]/(treesize*treesize);
	    for(k=0;k<4;k++){
	      q[k] /= treesize;
	      if(fabs(q[k]-eq[k]) > dbksran[l]) 
		dbksran[l] = fabs(q[k]-eq[k]);
	    }
	  }
	delete [] randsel;
	randsel=NULL;
      }

      for(l=0;l<iter;l++)
	if(dbksran[l]>dbks[i]) pvalue[i]+=1.0;
      pvalue[i] /= iter;

      delete [] seedran;
      seedran=NULL;
      delete [] seed;
      seed=NULL;
      delete [] tree;
      tree=NULL;
  }
  
  output <<fixed<<setprecision(8);
  output <<setw(7)<<"label"<<setw(12)<<"dbks"<<setw(12)<<"pvalue"<<endl;
  for(i=0;i<unique;i++){
    output <<setw(7)<<spuni[i];
    if(spused[i] == 0)
      output <<setw(12)<<"NA"<<setw(12)<<"NA"<<endl;
    else
      output <<setw(12)<<dbks[i]<<setw(12)<<pvalue[i]<<endl;
  }
  output.close();

  return 0;
}

// Helper function to the sort function. Sorts down an individual entry. From Numerical Recipes.  
void sift_down(sample tosort[], const int l, const int r)
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

// Sorting function, from Numerical Recipes.
void hpsort(sample tosort[])
{
  int i,n;
  n = sizeof(tosort)/sizeof(*tosort);
  for (i=n/2-1; i>=0; i--) sift_down(tosort, i, n-1);
  for (i=n-1; i>0; i--){
    swap(tosort[0],tosort[i]);
    sift_down(tosort,0,i-1);
  }
}

// Function that will generate an array of integers randomly ordered and containing the numbers
// from 0 to arraysize-1.
int *randgen(int length)
{
  int *arr;
  int randnum,j;
  arr = new int[length];
  arr[0] = 0;
  for(j=1;j<length;j++){
    randnum = (int) runif(0,j); 
    arr[j] = arr[randnum];
    arr[randnum] = j;
  }

  return arr; 
  
}

/* A function that uses ran1 to generate a random number between min and max (inclusive).
   Can be called as many times as ran1 will allow. Probably something like 100e6.*/      
float runif(float min, float max)
{

  double rannum,interval;

  if(idum == 0){
    int temp = (unsigned int) time( NULL );
    idum = -temp;
  }

  if(idum == 0){
    cout << "idum == 0, not safe to call ran1!\n";
    exit(1);
  }else{
    rannum = ran1();
  }
  interval = max - min;
  interval += 6.0e-15;
  rannum = min + (rannum*interval)-3.0e-15;
  if(rannum < min){ return min;
  }else if(rannum > max){ return max;
  }else{ return rannum;
  }
}

/* Random number generator from Numerical Recipes in C++ page 284*/

double ran1()
{

  const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
  const int NDIV=(1+(IM-1)/NTAB);
  const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
  static int iy=0;
  static int iv[NTAB];
  int j,k;
  double temp;

  if(idum<=0||!iy){
    if(-idum<1) idum=1;
    else idum=-idum;
    for(j=NTAB+7;j>=0;j--){
      k=idum/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if(idum<0)idum += IM;
      if(j<NTAB)iv[j]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if(idum<0)idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=idum;
  if((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
}
