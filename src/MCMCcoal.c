/* MCMCcoal.c
   Markov chain Monte Carlo coalescent program for population genetics data.
   It currently implements the JC69 model for analysing sequence data from multiple 
   loci.
   The input includes the number of species, the species tree, species 
   divergence times as well as theta for current and ancestral populations.
   The program then simulates a genealogical tree and evolves sequences along 
   it under the JC69 model.  Different gene trees are thus used for different 
   replicates.

   Copyright by Ziheng Yang, since July 2002

   cc -o MCMCcoal -O3 MCMCcoal.c tools.c -lm
   cc -o MCMCcoal  -m64 -march=opteron -mtune=opteron -ansi -O3 -funroll-loops -fomit-frame-pointer -finline-functions MCMCcoal.c tools.c -lm
   cc -o MCMCcoal -march=athlon -mcpu=athlon -O4 -funroll-loops -fomit-frame-pointer -finline-functions MCMCcoal.c tools.c -lm

   cc -o MCMCcoal -mcpu=G5 -O4 -funroll-loops -fomit-frame-pointer -finline-functions MCMCcoal.c tools.c -lm

   cl -O2 MCMCcoal.c tools.c
*/

#ifndef _MCMCCOAL_C_
#define _MCMCCOAL_C_

#include "MCMCcoal.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#define FPN(file) fputc('\n', file)
#define NSPECIES      200       /* max # of species */
#define NS            1000      /* max # of sequences */
#define NBRANCH       NS*2-2
#define MAXNSONS      2
#define NGENE         50000     /* max # of loci */
#define LSPNAME       50        /* # characters in sequence names */
//From paml.h
#define square(a) ((a)*(a))
#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define F0 stdout
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define Pi  3.1415926535897932384626433832795
#define DGammaMean 1
#define rndexp(mean) (-(mean)*log(rndu()))
#define FAST_RANDOM_NUMBER
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

struct CommonInfo {
   char *z[2*NS-1], *spname[NS], seqf[256], locusratef[256], cleandata, oldconP[NS*2-1];
   int ns, ls, lgene[1], model, clock, mcmc;
   int seqtype, ngene, *pose, npatt, np, readpattern;
   int ntime, ncatG, ncode, print, fix_rgene, posG[1];
   double alpha, pi[4], piG[1][4], *rates, rgene[1];
   double *conP;  /* not used */
   int    sconP, curconP;
   double *conPin[2], *conP0, *fpatt, space[1000000];
}  com;

struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[2], ibranch, ipop;  /* ibranch not used? */
   double branch, age, *conP, label;       /* age is uptodate but not branch */
   char fix_age;                       /* not used */
}  *nodes, *gnodes[NGENE], nodes_t[2*NS-1]; 
/*  nodes is a pointer, gnodes holds the gene trees, nodes_t is temp space
    gnodes is not allocated for the simulation program 
*/

/* sptree.nseqsp is working space, real thing in data.nseqsp[locus].
   npop is the number of theta's.  
   pops[] holds the node numbers in the species tree of the npop populations;
   the ancestral populations should be after the recent populations.  
   pop_pop_table[i][j] = 1 if pop i can coalescence in pop j and 0 if otherwise.
*/
struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nseqsp[NSPECIES]; 
   int npop, pops[NSPECIES*2-1], pop_pop_table[NSPECIES*2-1][NSPECIES*2-1];
   struct TREESPN {
      char name[LSPNAME*2];
      int father, nson, sons[2];
      double age, theta;
   } nodes[2*NSPECIES-1];
}  sptree;

enum {PrBranch=1, PrNodeNum=2, PrLabel=4, PrAge=8} OutTreeOptions;
enum {BASEseq=0, CODONseq, AAseq, CODON2AAseq} DataTypes;

int ReadSpeciesTree(char *datafile);
int DownSptreeSetPops(int inode);
void SimulateData(void);
void GetRandomGtree(int locus);
int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf,double space[]);
int abyx (double a, double x[], int n)
{ int i; for (i=0; i<n; x[i]*=a,i++) ;  return(0); }
void EvolveJC (int inode);
//void PMismatch3s (void);
void p0124Fromb0b1 (int itree, double p[5], double b[2]);
int ResetSpeciesGeneTree(int locus);
int Coalescence1Pop(int ispecies);
FILE *gfopen(char *filename, char *mode);
void starttime (void);
static void MCMCSetSeed (unsigned int seed);
void error2 (char * message);
int ReadaTreeN (FILE *ftree, int *haslength, int *haslabel, int copyname, int popline);
int OutaTreeN (FILE *fout, int spnames, int printopt);
static double rndu (void);
void printSeqs(FILE *fout, int *pose, char keep[], int format);
int MultiNomial2 (int n, int ncat, double prob[], int nobs[], double space[]);
int rndpoisson (double m);
static double LnGamma (double x);
int xtoy (double x[], double y[], int n);
int OutSubTreeN (FILE *fout, int inode, int spnames, int printopt, char *labelfmt);
static int DiscreteGamma (double freqK[], double rK[], double alpha, double beta, int K, int mean);
int MultiNomialAlias (int n, int ncat, double F[], int L[], int nobs[]);
int MultiNomialAliasSetTable (int ncat, double prob[], double F[], int L[], double space[]);
static double rndgamma (double s);
void ClearNode (int inode);
int ReadaTreeB (FILE *ftree, int popline);
int IsNameNumber(char line[]);
int print1seq (FILE*fout, char *z, int ls, int encoded, int pose[]);
long factorial(int n);
static double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
int matIout (FILE *fout, int x[], int n, int m);
void BranchToNode (void);
char *getcodon (char codon[], int icodon);
static double PointChi2 (double prob, double v);
static double PointNormal (double prob);

int noisy;
char timestr[32];
static double OLDAGE=99;
static int LASTROUND=1, debug=0, testlnL=0, totalTimes, totalSeqs, NPMat=0;
static char BASEs[]="TCAGUYRMKSWHBVDN?-";
static char AAs[] = "ARNDCQEGHILKMFPSTWYV*-?x";
static char BINs[] ="TC";

/*
int GetMem(void);
void FreeMem(void);
int ReadSeqData(char *seqfile, char *locusratef, FILE*fout, char cleandata);
double lnpGB_thetatau(int locus);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
int UseLocus(int locus, int copytreeconP, int useData, int setSeqName);
int AcceptLocus(int locus, int copyconP);
int GetInitials(void);
int collectx(double x[]);
int MCMC(FILE* fout);
int UpdateGB_InternalNode(double* lnL, double finetune);
int UpdateGB_SPR(double* lnL, double finetune);
int UpdateTheta(double finetune, double space[]);
int UpdateTau(double *lnL, double finetune, double space[]);
int mixing(double* lnL, double finetune, double space[]);
int UpdateLocusrateHeredity(double* lnL, double finetune);
int UpdateSequencingErrors(double* lnL, double finetune, double space[]);
void GraftNode(int source, int target, double age, int ipop);
int MatchGTree(void);
void printGtree(int printBlength);
void checkGtree(void);
int getab_beta(void);
*/



#define REALSEQUENCE
#define NODESTRUCTURE
/*#include "treesub.c"*/

void mcmccoal(char* simul_file_name)
{
   char datafile[128];
   char outfile[128]="out.txt", VStr[32]="Version 1.2, March 2007";
   //FILE *fout=(FILE*)gfopen(outfile,"w");
   int i;
   noisy=0;
   //printf("MCMCcoal %s\n", VStr);
   starttime();
   strcpy(datafile, simul_file_name); // [SJ] Changed
   com.cleandata=0;
   com.clock=1; com.ncode=4; com.model=0;
   FOR(i,4) com.pi[i]=1./4;
   ReadSpeciesTree(datafile);
   com.mcmc=0;
   SimulateData();
	 //fclose(fout);
}


int ReadSpeciesTree (char* datafile)
{
/* this reads the control data file for simulation.
   It reads the species tree and initializes sptree.
   Note that (sptree.nseqsp[j] > 1) adds a population.
   This also reads the gamma prior parameters data.a_gamma[] and data.b_gamma[].
*/
   FILE *fctl=(FILE*)gfopen(datafile,"r");
   char line[16000]={'\0'}, *pline, *seqerrstr="0EF";
   int i,j, lline=16000, haserr;
   double a,b, percent[2];
   //FILE*fseed=(FILE*)gfopen("SeedUsed", "w");

   fscanf(fctl, "%s%d", com.seqf, &i);
   if(i<0) i=abs(2*(int)time(NULL)+1);
   //fprintf(fseed, "%d\n", i); fclose(fseed);

   MCMCSetSeed(i);
   fscanf(fctl, "%d", &com.ns);
   if((sptree.nspecies=com.ns)>NSPECIES) error2("raise NSPECIES");

   nodes=nodes_t;  /* species tree is read into the temp space nodes_t */
   for(i=0; i<NS; i++) {  /* for both species & gene trees */
      if (com.spname[i]) free(com.spname[i]);
      com.spname[i]=(char*)malloc((LSPNAME+1)*sizeof(char));
      /* FOR(j,LSPNAME) com.spname[i][j]=0; */
   }
   for(i=0; i<sptree.nspecies; i++) {
      fscanf(fctl,"%s", com.spname[i]);
      strcpy(sptree.nodes[i].name, com.spname[i]);
   }
   
   for(i=0,j=0; i<sptree.nspecies; i++) {
      fscanf(fctl,"%d",&sptree.nseqsp[i]);
      j += sptree.nseqsp[i];
      if(sptree.nseqsp[i]<1) 
         puts("need at least 1 seq from each species?");
   }
   if(j>NS) { printf("\n%d sequences, too many.  Raise NS?", j);  exit(-1); }
   fgets(line, lline, fctl);
   if(sptree.nspecies>1) {
      ReadaTreeN (fctl, &i, &j, 0, 1);
      //OutaTreeN(F0,1,1);  FPN(F0);  //Do not print the tree [SJ]
   }
   else {
      tree.root=nodes[0].nson=0;  tree.nnode=1; tree.nbranch=0; 
      nodes[0].age=0;
   }
/*
   printf("%d species: ", sptree.nspecies);
   for(i=0; i<sptree.nspecies; i++) 
      printf(" %s (%d)",com.spname[i],sptree.nseqsp[i]);
*/
   /* copy into sptree */
   sptree.nnode=tree.nnode;  sptree.nbranch=tree.nbranch; sptree.root=tree.root;
   for(i=0; i<sptree.nspecies*2-1; i++) {
      sptree.nodes[i].father = nodes[i].father;
      sptree.nodes[i].nson = nodes[i].nson;
      for(j=0;j<sptree.nodes[i].nson;j++) 
         sptree.nodes[i].sons[j] = nodes[i].sons[j];
      sptree.nodes[i].age = nodes[i].age = nodes[i].branch;
      sptree.nodes[i].theta = nodes[i].label;
   }

   for(i=0,sptree.npop=0; i<sptree.nspecies; i++)
      if(sptree.nseqsp[i]>1)  sptree.pops[sptree.npop++]=i;
   for(; i<2*sptree.nspecies-1; i++) sptree.nodes[i].name[0]='\0';

   DownSptreeSetPops(sptree.root);

//   printf("\npopulation by population table:\n");
   for(i=0; i<2*sptree.nspecies-1; i++) for(j=0; j<2*sptree.nspecies-1; j++) 
      sptree.pop_pop_table[i][j]=0;
   for(i=0; i<2*sptree.nspecies-1; i++) {
      for(j=i; ; j=sptree.nodes[j].father) {
         if(j>=sptree.nspecies || sptree.nseqsp[j]>1)
            sptree.pop_pop_table[i][j]=1;
         if(j==sptree.root) break;
      }
   }
/*   for(i=0; i<2*sptree.nspecies-1; i++,FPN(F0)) {
      printf("species %2d %-10s ", i+1,sptree.nodes[i].name);
      printf("age %.5f  theta %.5f  ", sptree.nodes[i].age, sptree.nodes[i].theta);
      for(j=0; j<2*sptree.nspecies-1; j++) 
         printf("%1d", sptree.pop_pop_table[i][j]);
      if(i<sptree.nspecies && sptree.nodes[i].age) 
         printf("\a  <-- age>0??");
   }
*/
   com.np = sptree.npop+sptree.nspecies-1;
/*
   if(!com.mcmc) {
      printf("\n%d theta parameters (populations) in the order:", sptree.npop);
      for(i=0; i<sptree.npop; i++)
         printf(" %2d (%s)", sptree.pops[i]+1, sptree.nodes[sptree.pops[i]].name);
      if(sptree.nspecies>1) {
         printf("\n%d species divergence times in the order:", sptree.nspecies-1);
         for(i=sptree.nspecies; i<sptree.nspecies*2-1; i++)
            printf(" %2d (%s)", i+1, sptree.nodes[i].name);
      }
   }
*/
   /* count the number of sequences in gene tree com.ns or data.ns[] */
   for(i=0,com.ns=0; i<sptree.nspecies; i++) com.ns+=sptree.nseqsp[i];
   if(com.ns>NS) error2("raise NS in source");
 	 fclose(fctl);
  return(0);
}

int DownSptreeSetPops (int inode)
{
/* This traverse down the species tree to see which nodes are visited,
   and initializes npop, pops[], and sptree.nodes[].name.
   Populations where coalescent is possible (theta's) are collected into 
   sptree.pops[] in the post-order tree traversal.
*/
   int k,ison;

   for (k=0; k<sptree.nodes[inode].nson; k++) {
      ison=sptree.nodes[inode].sons[k];
      if(sptree.nodes[ison].nson)  
         DownSptreeSetPops(ison);
      if(sptree.nspecies<11)   /* danger of memory overrun!!! */
         strcat(sptree.nodes[inode].name, sptree.nodes[ison].name);
   }
   if(inode>=sptree.nspecies)
      sptree.pops[sptree.npop++]=inode;
   return(0);
}

void SimulateData (void)
{
   char *treef="out.trees", savetree=1,makeseq=0; // makeseq=0 means that the program does not simulate sequences [SJ] 
   //FILE *fseq=gfopen(com.seqf, "w");
	 FILE	*ftree=gfopen(treef,"w");
   int nr, ir, i,j, variable_ns=0, nseqsp0[NSPECIES];
   double Tree1=0, mtH3=0, mtmrca=0;

   /* PMismatch3s(); */
   nr=1; com.ls=1000; com.cleandata=1;
//   printf("\nnloci & ls? ");
//   scanf("%d%d", &nr, &com.ls);

   com.alpha=0;  com.ncatG=8; for(i=0;i<4;i++) com.pi[i]=1./4;
//   printf("\nsequence data go into %s, trees into %s.\n", com.seqf, treef);
//   com.z[0]=(char*)malloc((2*com.ns-1)*com.ls*sizeof(char));
//   for (i=1; i<2*com.ns-1; i++) com.z[i]=com.z[i-1]+com.ls;
//   if(com.alpha) com.rates=(double*)malloc(com.ls*sizeof(double));
//   if(com.z[0]==NULL || (com.alpha && com.rates==NULL)) { printf("\nError: oom.\n"); exit(-1); };

   if(variable_ns) FOR(j,sptree.nspecies) nseqsp0[j]=sptree.nseqsp[j];
   //if(makeseq) fprintf(fseq, "%6d\n", nr);
   for (ir=0; ir<nr; ir++) {
      if(variable_ns)
         for( ; ; ) { /* make sure at least 2 sequences are in the data */
            for(j=0,com.ns=0;j<sptree.nspecies;j++) 
               com.ns+=sptree.nseqsp[j]=(int)(nseqsp0[j]*rndu()*0.95+0.5);
            if(com.ns>1) break;
         }
      GetRandomGtree(-1);
      if(savetree) { OutaTreeN(ftree,1,1);  FPN(ftree); }
/*
if(MatchGTree()) {
   Tree1++;
   mtH3=(mtH3*(Tree1-1)+nodes[nodes[2].father].age)/Tree1;
   mtmrca=(mtmrca*(Tree1-1)+nodes[tree.root].age)/Tree1;
}
if((ir+1)%10000==0) 
   printf("\r%d  %9.6f\tEh3%9.6f EtMRCA %9.6f", ir+1, Tree1/(ir+1),mtH3,mtmrca);
*/
/*      if(!makeseq) continue;
      if (com.alpha) {
         Rates4Sites (com.rates, com.alpha, com.ncatG, com.ls, 1, com.space);
         for (j=1; j<com.ls; j++) com.rates[j]+=com.rates[j-1];
         abyx (1/com.rates[com.ls-1], com.rates, com.ls);
      }
      FOR(i,com.ls) com.z[tree.root][i]=(char)(rndu()*4.);
      EvolveJC (tree.root);
*/
//      for(i=0;i<sptree.nspecies;i++) /*fprintf(fseq," %3d",sptree.nseqsp[i])*/; 
      //FPN(fseq);
      //printSeqs(fseq, NULL, NULL, 0);

   }  /* for(ir) */
   //fclose(fseq);
//	 free(com.z[0]);
   if(com.alpha) free(com.rates);
   if(variable_ns) FOR(j,sptree.nspecies) sptree.nseqsp[j]=nseqsp0[j];
	 fclose(ftree);
   //exit(0);
}

/* used by Coalescence1Pop(). */
static int cNodeNumber, nin_G[NSPECIES*2-1], ins_G[NSPECIES*2-1][NS];

void GetRandomGtree(int locus)
{
/* generates a random gene tree in nodes[].  This does not initialize ipop.
*/
   int i;

   ResetSpeciesGeneTree(locus);
   Coalescence1Pop(sptree.root);
   /* NodeToBranch(); */
   tree.root=cNodeNumber-1;  nodes[tree.root].branch=0;  
   nodes[tree.root].father=-1;  tree.nnode=com.ns*2-1;
   FOR (i,tree.nnode) if(i!=tree.root) 
      nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
}

static int n123marks[][3]= {{1,2,3},{1,2,3},{2,3,1},{3,1,2}};

void p0124Fromb0b1 (int itree, double p[5], double b[2])
{
/* This calculates p0,p1,p2,p3,p4 for the 5 site patterns for 3 species, 
   given branch lengths b0 and b1.  b0 is the gap, and b1 is the distance 
   from the ancestor of 1 and 2 to species 1.
*/
   double e1, e2, e3;

   e1=exp(-4./3*b[1]);  e2=exp(-8./3*(b[0]+b[1]));  e3=e1*e2;  e1=e1*e1;
   p[0]      = (1. + 3*e1 +  6*e2 +  6*e3)/16;
   p[n123marks[itree][0]]  = (3. + 9*e1 -  6*e2 -  6*e3)/16;
   p[n123marks[itree][1]]=p[n123marks[itree][2]] 
      = (3. - 3*e1 +  6*e2 -  6*e3)/16;
   p[4]      = (6. - 6*e1 - 12*e2 + 12*e3)/16;
}


void EvolveJC (int inode)
{
/* Special version of Evolve that works with JC69-like (poisson) model only.
   This can be used to generate amino acid sequences also.
   For each branch in the tree, this determines the number of mutations 
   and then assigns mutations to sites at random.
   When alpha>0, com.rates[] are the accumulative probabilities.
*/
   int is,j, h, nmut, imut, ison;
   double r;

   if (com.alpha && fabs(com.rates[com.ls-1]-1)>1e-4) 
      { printf ("rates c.d.f.: 1 = %.6f?\n", com.rates[com.ls-1]); exit(-1); }
   for (is=0; is<nodes[inode].nson; is++) {
      ison=nodes[inode].sons[is];
      FOR(h,com.ls) com.z[ison][h]=com.z[inode][h];
      nmut = rndpoisson (nodes[ison].branch*com.ls);
      for (imut=0; imut<nmut; imut++) {
         if (com.alpha==0) h=(int)(rndu()*com.ls);
         else 
            for (h=0,r=rndu(); h<com.ls; h++) if (r<com.rates[h]) break;
         j=(int)(rndu()*3);
         if (j>=com.z[ison][h]) j++;
         com.z[ison][h]=(char)j;
      }
      if (nodes[ison].nson) EvolveJC(ison);
   }
}

int OutaTreeN (FILE *fout, int spnames, int printopt)
{
/* 
   print the current tree.
*/
   int i, intlabel=1;
   char* labelfmt[2]={"#%.5f", "#%.0f"};

   if(printopt & PrLabel) {
      for(i=0;i<tree.nnode;i++) 
         if(nodes[i].label-(int)nodes[i].label != 0) intlabel=0;
   }

   OutSubTreeN(fout,tree.root,spnames,printopt, labelfmt[intlabel]);
   if(printopt&PrNodeNum) fprintf(fout," %d ", tree.root+1);
   if((printopt & PrLabel) && nodes[tree.root].label>0) 
      fprintf(fout, labelfmt[intlabel], nodes[tree.root].label);
   if((printopt & PrAge)  && nodes[tree.root].age) 
      fprintf(fout, " @%.3f", nodes[tree.root].age);

   if((printopt&PrBranch) && nodes[tree.root].branch>1e-8)
      fprintf(fout,": %.6f", nodes[tree.root].branch);

   fputc(';',fout);
	 fflush(fout);
   return(0);
}

int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf,
    double space[])
{
/* Rates for sites from the gamma (ncatG=0) or discrete-gamma (ncatG>1).
   Rates are converted into the c.d.f. if cdf=1, which is useful for
   simulation under JC69-like models. 
   space[ncatG*5]
*/
   int h, ir, j, K=ncatG, *Lalias=(int*)(space+3*K), *counts=(int*)(space+4*K);
   double *rK=space, *freqK=space+K, *Falias=space+2*K;

   if (alpha==0) 
      { if(rates) FOR(h,ls) rates[h]=1; }
   else {
      if (K>1) {
         DiscreteGamma (freqK, rK, alpha, alpha, K, DGammaMean);

         MultiNomialAliasSetTable(K, freqK, Falias, Lalias, space+5*K);
         MultiNomialAlias(ls, K, Falias, Lalias, counts);

         for (ir=0,h=0; ir<K; ir++) 
            for (j=0; j<counts[ir]; j++)  rates[h++]=rK[ir];
      }
      else 
         for (h=0; h<ls; h++) rates[h]=rndgamma(alpha)/alpha;
      if (cdf) {
         for (h=1; h<ls; h++) rates[h]+=rates[h-1];
         abyx (1/rates[ls-1], rates, ls);
      }
   }
   return (0);
}


int ResetSpeciesGeneTree (int locus)
{
/* Resets the species tree for Coalescence1Pop(), initializing nin_G, 
   ins_G[], ipop.
   Also resets cNodeNumber for constructing the gene tree.
   nodes[].ipop for tips in the gene tree are set already.
*/
   int is,j;

#ifdef MetropolisHastings
   com.ns=data.ns[locus];
   nodes=gnodes[locus]; 
   for(is=0; is<sptree.nspecies; is++) sptree.nseqsp[is]=data.nseqsp[locus][is];
#endif
   for(is=0; is<sptree.nspecies; is++) 
      nin_G[is] = sptree.nseqsp[is];
   for(is=0,cNodeNumber=0; is<sptree.nspecies; is++) {
      for(j=0; j<sptree.nseqsp[is]; j++,cNodeNumber++) {
         ins_G[is][j]=cNodeNumber;
// [SJ] Changes in version 1.03
//         if(nin_G[is]>1) 
            sprintf(com.spname[cNodeNumber], "%s%d", sptree.nodes[is].name,j+1);
//         else 
//            sprintf(com.spname[cNodeNumber], "%s", sptree.nodes[is].name);
      }
   }
   for(; is<2*sptree.nspecies-1; is++) nin_G[is]=0;
   if(cNodeNumber!=com.ns)  error2("ns mismatch");
   for(j=0;j<com.ns;j++) ClearNode(j);
   return(0);
}

int Coalescence1Pop (int ispecies)
{
/* This simulates the coalescent process in population or species ispecies.
   It generates the random genealogy tree (possibly consisting of several 
   disconnected subtrees) with waiting times.
   t: waiting time; T: node age
   This is used by GetRandomGtree().  
   nin_G[] and ins_G[] are set in ResetSpeciesGeneTree().
*/
   int j, k,k1,k2, father=sptree.nodes[ispecies].father;
   double t, T;

   for (k=0; k<sptree.nodes[ispecies].nson; k++)
      Coalescence1Pop(sptree.nodes[ispecies].sons[k]);
   T=sptree.nodes[ispecies].age;
   if(nin_G[ispecies]>1 && sptree.nodes[ispecies].theta==0)
      error2("theta=0 for a population");
   for (j=nin_G[ispecies]; j>1; j--,cNodeNumber++) {
      t=rndexp(sptree.nodes[ispecies].theta/(j*(j-1.)));  T += t;
      if(ispecies!=sptree.root && T>sptree.nodes[father].age) break;

      k  = (int)(j*rndu());  
      k1 = ins_G[ispecies][k]; ins_G[ispecies][k] = ins_G[ispecies][j-1];
      k  = (int)((j-1)*rndu());
      k2 = ins_G[ispecies][k]; ins_G[ispecies][k] = cNodeNumber;

      nodes[cNodeNumber].nson = 2;
      nodes[cNodeNumber].ipop = ispecies;
      nodes[cNodeNumber].sons[0] = k1;
      nodes[cNodeNumber].sons[1] = k2;
      nodes[k1].father=nodes[k2].father = cNodeNumber;
      nodes[cNodeNumber].age = T;
   }
  /* outgoing lineages added to father's list */
   nin_G[ispecies] = j;
   if(ispecies==sptree.root && j>1) 
      error2("nin_Gtree > 1 at root");
   if(ispecies!=sptree.root) {
      for(k=0; k<j; k++) 
         ins_G[father][nin_G[father]++] = ins_G[ispecies][k];
   }
   return (0);
}


void printGtree (int printBlength)
{
   int i,j, ipop, ipopTrue;
   double t, tb[2];

   for(i=0; i<tree.nnode; i++) 
      if(i!=tree.root) {
         nodes[i].branch = nodes[nodes[i].father].age-nodes[i].age;
         if(nodes[i].branch<0) 
            printf("");
      }
   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%7s%7s%11s (%s) %7s%7s","father","node","time","ipop","nson:","sons");
   FOR (i, tree.nnode) {
      t=nodes[i].age;
      ipop=nodes[i].ipop;
      tb[0]=sptree.nodes[ipop].age;  tb[1]=OLDAGE;
      if(ipop!=sptree.root) tb[1]=sptree.nodes[sptree.nodes[ipop].father].age;

      ipopTrue=(t>=tb[0] && t<tb[1]);
      printf ("\n%7d%7d %11.6f (%2d %s %11.6f):%7d  ",
         nodes[i].father, i, t, ipop, (ipopTrue?"true ":"false"), 
         sptree.nodes[ipop].age, nodes[i].nson);
      FOR(j, nodes[i].nson) printf (" %2d", nodes[i].sons[j]);
   }
   FPN(F0); OutaTreeN(F0,0,0); FPN(F0); OutaTreeN(F0,1,0); 
   if(printBlength) {
      FPN(F0); OutaTreeN(F0,1,1); FPN(F0); 
   }
}


int MatchGTree (void)
{
/* This tests the gene tree topology
*/
   return 
     (nodes[0].father==nodes[1].father
   && nodes[3].father==nodes[4].father 
   && nodes[nodes[3].father].father==tree.root
   && nodes[nodes[0].father].age<sptree.nodes[sptree.root].age
   && nodes[nodes[2].father].age<sptree.nodes[sptree.root].age
   && nodes[nodes[3].father].age<sptree.nodes[sptree.root].age
   );  /* P(gtree) = 0.237913 */
}

FILE *gfopen(char *filename, char *mode)
{
   FILE *fp;

   if(filename==NULL || filename[0]==0) 
      error2("file name empty.");

   fp=(FILE*)fopen(filename, mode);
   if(fp==NULL) {
      printf("\nerror when opening file %s\n", filename);
      if(!strchr(mode,'r')) exit(-1);
      printf("tell me the full path-name of the file? ");
      scanf("%s", filename);
      if((fp=(FILE*)fopen(filename, mode))!=NULL)  return(fp);
      puts("Can't find the file.  I give up.");
      exit(-1);
   }
   return(fp);
}

static time_t time_start;

void starttime (void)
{
   time_start=time(NULL);
}

static unsigned int z_rndu=137;
static unsigned int w_rndu=123456757;

void MCMCSetSeed (unsigned int seed)
{
   if(sizeof(int) != 4) 
      puts("oh-oh, we are in trouble.  int not 32-bit?");
   z_rndu=170*(seed%178)+137;
   w_rndu = seed*127773;
}

void error2 (char * message)
{ printf("\nError: %s.\n", message); exit(-1); }

static int *CladeLabel=NULL;

int PopPaupTreeRubbish(FILE *ftree);

int PopPaupTreeRubbish(FILE *ftree)
{
/* This reads out the string in front of the tree in the nexus format, 
   typically "tree PAUP_1 = [&U]" with "[&U]" optional
*/
   int ch;

   for (; ;) {
      ch=fgetc(ftree);
      if(ch=='(')          { ungetc(ch,ftree); return(0); }
      else if(ch==EOF)     return(-1);
   }
   return(0);
}

int ReadaTreeN (FILE *ftree, int *haslength, int *haslabel, int copyname, int popline)
{
/* Read a tree from ftree, using the parenthesis node representation of trees.
   Branch lengths are read in nodes[].branch, and branch (node) labels 
   (integers) are preceeded by # and read in nodes[].label.  If the clade label
   $ is used, the label is read into CladeLabel[] first and then moved into
   nodes[].label in the routine DownTreeCladeLabel().

   This assumes that com.ns is known.
   Species names are considered case-sensitive, with trailing spaces ignored.

   copyname = 0: species numbers and names are both accepted, but names have 
                 to match the names in com.spname[], which are from the 
                 sequence data file.  Used by baseml and codeml, for example.
              1: species names are copied into com.spname[], but species 
                 numbers are accepted.  Used by evolver for simulation, 
                 in which case no species names were read before.
              2: the tree must have species names, which are copied into com.spname[].
                 Note that com.ns is assumed known.  To remove this restrition, 
                 one has to consider the space for nodes[], CladeLabel, starting 
                 node number etc.

   isname = 0:   species number; 1: species name; 2:both number and name
*/
   int cnode, cfather=-1;  /* current node and father */
   int inodeb=0;  /* node number that will have the next branch length */
   int i,j, level=0, isname, ch=' ', icurspecies=0;
   char check[NS], delimiters[]="(),:#$=@><;", skips[]="\"\'";
   int lline=32000;
   char line[32000];

   if(com.ns<=0)  error2("need to know ns before reading tree.");

   if((CladeLabel=(int*)malloc((com.ns*2-1)*sizeof(int)))==NULL) 
      error2("Out of Memory... trying to get space for cladelabel");
   FOR(i,2*com.ns-1) CladeLabel[i]=0;

   /* initialization */
   FOR(i,com.ns) check[i]=0;
   *haslength=0, *haslabel=0;
   tree.nnode=com.ns;  tree.nbranch=0;
   FOR(i,2*com.ns-1) {
      nodes[i].father=nodes[i].ibranch=-1;
      nodes[i].nson=0;  nodes[i].label=0;  nodes[i].branch=0;
      nodes[i].age=0;  /* TipDate models set this for each tree later. */
#if (defined(BASEML) || defined(CODEML))
      nodes[i].fossil=0;
#endif
   }
   while(isspace(ch)) ch=fgetc(ftree);  /* skip spaces */
   ungetc(ch,ftree);
   if (isdigit(ch)) { ReadaTreeB(ftree,popline); return(0); }

   PopPaupTreeRubbish(ftree);

   for ( ; ; ) {
      ch = fgetc (ftree);
      if (ch==EOF) return(-1);
      else if (!isgraph(ch) || ch==skips[0] || ch==skips[1]) continue;
      else if (ch=='(') {       /* left (  */
         level++;
         cnode=tree.nnode++;
         if(tree.nnode>2*com.ns-1) 
			 error2("check #seqs and tree: perhaps too many '('?");
         if (cfather>=0) {
            if(nodes[cfather].nson>=MAXNSONS)
               error2("too many daughter nodes, raise MAXNSONS");
            nodes[cfather].sons[nodes[cfather].nson++] = cnode;
            nodes[cnode].father=cfather;
            tree.branches[tree.nbranch][0]=cfather;
            tree.branches[tree.nbranch][1]=cnode;
            nodes[cnode].ibranch=tree.nbranch++;
         }
         else
            tree.root=cnode;
         cfather=cnode;
      }
      /* treating : and > in the same way is risky. */
      else if (ch==')') { level--;  inodeb=cfather; cfather=nodes[cfather].father; }
      else if (ch==':'||ch=='>') { 
         if(ch==':') *haslength=1;
         fscanf(ftree,"%lf",&nodes[inodeb].branch); 
      }
      else if (ch=='#'||ch=='<') { *haslabel=1; fscanf(ftree,"%lf",&nodes[inodeb].label); }
      else if (ch=='$')          { *haslabel=1; fscanf(ftree,"%d",&CladeLabel[inodeb]); }
      else if (ch=='@'||ch=='=') { 
         fscanf(ftree,"%lf",&nodes[inodeb].age);
#if (defined(BASEML) || defined(CODEML))
         if(com.clock) nodes[inodeb].fossil=1;
#endif
      }
      else if (ch==',') ;
      else if (ch==';' && level!=0) error2("; in treefile");
      else { /* read species name or number */
         line[0]=(char)ch;  line[1]=(char)fgetc(ftree);
/*         if(line[1]==(char)EOF) error2("eof in tree file"); */

         for (i=1; i<lline; )  { /* read species name into line[] until delimiter */
            if ((strchr(delimiters,line[i]) && line[i]!='@') 
               || line[i]==(char)EOF || line[i]=='\n')
               { ungetc(line[i],ftree); line[i]=0; break; }
            line[++i]=(char)fgetc(ftree);
         }
         for(j=i-1;j>0;j--) /* trim spaces*/
            if(isgraph(line[j])) break; else line[j]=0;
         isname=IsNameNumber(line);

         if (isname==2) {       /* both number and name */
            sscanf(line, "%d", &cnode);   cnode--;
            strcpy(com.spname[cnode], line);
         }
         else if (isname==0) {  /* number */
            if(copyname==2) error2("Use names in tree.");
            sscanf(line, "%d", &cnode);   cnode--;
         }
         else {                 /* name */
            if(!copyname) {
               for(i=0; i<com.ns; i++) if (!strcmp(line,com.spname[i])) break;
               if((cnode=i)==com.ns) { printf("\nSpecies %s?\n", line); exit(-1); }
            }
            else {
               if(icurspecies>com.ns-1) {
                  error2("error in tree: too many species in tree");
               }
               strcpy(com.spname[cnode=icurspecies++], line);
            }
         }

         nodes[cnode].father=cfather;
         if(nodes[cfather].nson>=MAXNSONS)
            error2("too many daughter nodes, raise MAXNSONS");

         nodes[cfather].sons[nodes[cfather].nson++]=cnode;
         tree.branches[tree.nbranch][0]=cfather;
         tree.branches[tree.nbranch][1]=cnode;
         nodes[cnode].ibranch=tree.nbranch++;
         check[inodeb=cnode]++;
      }
      if (level<=0) break;
   }
   /* read branch length and label for the root if any */
   /* treating : and > in the same way is risky. */
   for ( ; ; ) {
      for(; ;) {
         ch=fgetc(ftree);
         if(ch=='\n' || ch==EOF || (isgraph(ch) && ch!=skips[0] && ch!=skips[1]))
            break;
      }
      if (ch==':'||ch=='>')       fscanf(ftree, "%lf", &nodes[tree.root].branch);
      else if (ch=='#'||ch=='<') { *haslabel=1; fscanf(ftree,"%lf",&nodes[inodeb].label); }
      else if (ch=='$')           error2("why do you label the whole tree?");
      else if (ch=='@'||ch=='=') { 
         fscanf(ftree,"%lf",&nodes[inodeb].age);
#if (defined(BASEML) || defined(CODEML))
         if(com.clock) nodes[inodeb].fossil=1;
#endif
      }
      else if (ch==';')  break;
      else  /* anything unrecognised */
         { ungetc(ch,ftree);  break; }
   }

   if (popline) fgets (line, lline, ftree);
   FOR(i,com.ns) {
      if(check[i]>1) 
         { printf("\nSeq #%d occurs more than once in the tree\n",i+1); exit(-1); }
      else if(check[i]<1) 
         { printf("\nSeq #%d (%s) is missing in the tree\n",i+1,com.spname[i]); exit(-1); }
   }
   if(tree.nbranch>2*com.ns-2) { 
      printf("nbranch %d", tree.nbranch); puts("too many branches in tree?");
   }
   if (tree.nnode != tree.nbranch+1) {
      printf ("\nnnode%6d != nbranch%6d + 1\n", tree.nnode, tree.nbranch);
      exit(-1);
   }

/* check that it is o.k. to comment out this line
   com.ntime = com.clock ? (tree.nbranch+1)-com.ns+(tree.root<com.ns)
                         : tree.nbranch;
*/

#if(defined(BASEML) || defined(CODEML))
   /* check and convert clade labels $ */
   if(com.clock>1 || (com.seqtype==1 && com.model>=2)) {
      for(i=0,j=0; i<tree.nnode; i++) {
         if(i<com.ns && CladeLabel[i]) 
            error2("cannot have $ label for tips.");
         if(CladeLabel[i]) j++;
      }
      if(j) /* number of clade labels */
         DownTreeCladeLabel(tree.root);

      /*
      OutaTreeN(F0,1,PrBranch|PrNodeNum);  FPN(F0);
      */

      for(i=0,com.nbtype=0; i<tree.nnode; i++) { 
         if(i==tree.root) continue;
         j=(int)nodes[i].label;
         if(j+1>com.nbtype)  com.nbtype=j+1;
         if(j<0||j>tree.nbranch-1)  error2("branch label in the tree");
      }
      if (com.nbtype<=1)
         error2("need branch labels in the tree for the model.");
      else {
         printf("\n%d branch types are in tree. Stop if wrong.", com.nbtype);
      }

#if(defined(CODEML))
      if(com.seqtype==1 && com.NSsites && (com.model==2 || com.model==3) && com.nbtype!=2)
         error2("only two branch types are allowed for branch models.");
#endif

   }
#endif

   free(CladeLabel);
   return (0);
}


#ifdef FAST_RANDOM_NUMBER
double rndu (void)
{
/* From Ripley (1987, page 46). 32-bit integer assumed */
   w_rndu = w_rndu*69069+1;
   return ldexp((double)w_rndu, -32);
}


void rndu_vector (double r[], int n)
{
/* From Ripley (1987, page 46). 32-bit integer assumed */
   int i;

   for(i=0; i<n; i++) {
      w_rndu = w_rndu*69069+1;
      r[i] = ldexp((double)w_rndu, -32);
   }
}


#else 

double rndu (void)
{
/* U(0,1): AS 183: Appl. Stat. 31:188-190 
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
*/
   static unsigned int x_rndu=11, y_rndu=23;
   double r;

   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
/*
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
*/
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}

#endif


void printSeqs(FILE *fout, int *pose, char keep[], int format)
{
/* Print sequences into fout, using paml (format=0) or paup (format=1) formats.
   Use pose=NULL if called before site patterns are collapsed.  
   keep[] marks the sequences to be printed.  Use NULL for keep if all sequences 
   are to be printed.
   Sequences may (com.cleandata==1) and may not (com.cleandata=0) be coded.
   com.z[] has site patterns if pose!=NULL.
   This uses com.seqtype, and com.ls is the number of codons for codon seqs.
   See notes in print1seq()

   format=0,1: PAML sites or patterns
          2:   PAUP Nexus format.


   This is used by evolver.  Check and merge with printsma().

*/
   int h, j,ls1, n31=(com.seqtype==1?3:1), nskept=com.ns, wname=30;
   const char *dt=(com.seqtype==AAseq?"protein":"dna");

   if(!com.cleandata)  
      error2("check different routines for printing seqs");

   ls1=(format==1?com.npatt:com.ls);
   if(keep) FOR(j,com.ns) nskept -= !keep[j];
   if(format==0 || format==1) 
      fprintf(fout, "\n\n%6d %7d %s\n\n", nskept, ls1*n31, (format==1?" P":""));
   else if (format==2) {  /* NEXUS format */
      fprintf(fout,"\nbegin data;\n");
      fprintf(fout,"   dimensions ntax=%d nchar=%d;\n", nskept,ls1);
      fprintf(fout,"   format datatype=%s missing=? gap=-;\n   matrix\n",dt);
   }

   for(j=0; j<com.ns; j++,FPN(fout)) {
      if(keep && !keep[j]) continue;
      fprintf(fout,"%s%-*s  ", (format==2?"      ":""),wname,com.spname[j]);
//      fprintf(fout,"%s%s", com.spname[j], "  "); // [SJ] v.1.03
      print1seq(fout,com.z[j], (format==1?com.npatt:com.ls), com.cleandata, pose);
   }
   if(format==2) fprintf(fout, "   ;\nend;");
   else if (format==1) {
       for(h=0; h<com.npatt; h++) {
/*
         fprintf(fout," %12.8f", com.fpatt[h]/(double)com.ls);
*/
         fprintf(fout," %4.0f", com.fpatt[h]);

         if((h+1)%15==0) FPN(fout);
      }
   }
   fprintf(fout,"\n\n");
   fflush(fout);
}


int MultiNomial2 (int n, int ncat, double prob[], int nobs[], double space[])
{
/* sample n times from a mutinomial distribution M(ncat, prob[])
   prob[] is considered cumulative prob if (space==NULL)
   ncrude is the number or crude categories, and lcrude marks the
   starting category for each crude category.  These are used 
   to speed up the process when ncat is large.
*/
   int i, j, crude=(ncat>20), ncrude, lcrude[200];
   double r, *pcdf=(space==NULL?prob:space), small=1e-5;

   ncrude=max2(5,ncat/20); ncrude=min2(200,ncrude);
   FOR(i,ncat) nobs[i]=0;
   if (space) {
      xtoy(prob, pcdf, ncat);
      for(i=1; i<ncat; i++) pcdf[i]+=pcdf[i-1];
   }
   if (fabs(pcdf[ncat-1]-1) > small) 
      error2("sum P!=1 in MultiNomial2");
   if (crude) {
      for(j=1,lcrude[0]=i=0; j<ncrude; j++)  {
         while (pcdf[i]<(double)j/ncrude) i++;
         lcrude[j]=i-1;
      }
   }
   FOR(i,n) {
      r=rndu();
      j=0;
      if (crude) {
         for (; j<ncrude; j++) if (r<(j+1.)/ncrude) break;
         j=lcrude[j];
      }
      for (; j<ncat-1; j++) if (r<pcdf[j]) break;
      nobs[j] ++;
   }
   return (0);
}     


int rndpoisson (double m)
{
/* m is the rate parameter of the poisson
   Numerical Recipes in C, 2nd ed. pp. 293-295
*/
   static double sq, alm, g, oldm=-1;
   double em, t, y;

/* search from the origin
   if (m<5) { 
      if (m!=oldm) { oldm=m; g=exp(-m); }
      y=rndu();  sq=alm=g;
      for (em=0; ; ) {
         if (y<sq) break;
         sq+= (alm*=m/ ++em);
      }
   }
*/
   if (m<12) { 
      if (m!=oldm) { oldm=m; g=exp(-m); }
      em=-1; t=1;
      for (; ;) {
         em++; t*=rndu();
         if (t<=g) break;
      }
   }
   else {
     if (m!=oldm) {
        oldm=m;  sq=sqrt(2*m);  alm=log(m);
        g=m*alm-LnGamma(m+1);
     }
     do {
        do {
           y=tan(3.141592654*rndu());
           em=sq*y+m;
        } while (em<0);
        em=floor(em);
        t=0.9*(1+y*y)*exp(em*alm-LnGamma(em+1)-g);
     } while (rndu()>t);
   }
   return ((int) em);
}


double LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.

   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x-1;

   if((double)nx==x && nx>0 && nx<10)
      lng=log((double)factorial(nx));
   else {
      if(x<=0) {
         error2("lnGamma not implemented for x<0");
         if((int)x-x==0) { puts("lnGamma undefined"); return(-1); }
         for (fneg=1; x<0; x++) fneg/=x;
         if(fneg<0) error2("strange!! check lngamma");
         fneg=log(fneg);
      }
      if (x<7) {
         f=1;  z=x-1;
         while (++z<7)  f*=z;
         x=z;   f=-log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673 
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}

int xtoy (double x[], double y[], int n)
{ int i; for (i=0; i<n; y[i]=x[i],i++) ;  return(0); }

int OutSubTreeN (FILE *fout, int inode, int spnames, int printopt, char *labelfmt)
{
   int i,ison;

   fputc ('(', fout);
   for(i=0; i<nodes[inode].nson; i++) {
      ison=nodes[inode].sons[i];
      if(nodes[ison].nson==0) {
         if(spnames) {
            if(printopt&PrNodeNum) fprintf(fout, "%d_",ison+1);
            fprintf(fout,"%s",com.spname[ison]);
         }
         else 
            fprintf(fout,"%d",ison+1);
      }
      else
         OutSubTreeN(fout, ison, spnames, printopt, labelfmt);

      if((printopt & PrNodeNum) && nodes[ison].nson) 
         fprintf(fout," %d ",ison+1);   /* fprintf(fout,"%d",ison+1-sptree.nspecies);  */
      if((printopt & PrLabel) && nodes[ison].label>0)
         fprintf(fout, labelfmt, nodes[ison].label);
      if((printopt & PrAge) && nodes[ison].age) 
         fprintf(fout, " @%.3f", nodes[ison].age);


/*  Add branch labels to be read by Rod Page's TreeView. */
#if (defined CODEML)
      if(printopt & PrLabel)
         fprintf(fout," #%.4f ", nodes[ison].omega);
#elif (defined EVOLVER)
      if((printopt & PrLabel) && nodes[ison].nodeStr) 
         fprintf(fout," '%s'", nodes[ison].nodeStr);
#endif
      if(printopt & PrBranch) fprintf(fout,": %.6f", nodes[ison].branch);

      if(i<nodes[inode].nson-1) fprintf(fout,", ");
   }
   fputc (')', fout);
   return (0);
}


int DiscreteGamma (double freqK[], double rK[], double alpha, double beta, int K, int mean)
{
/* discretization of gamma distribution with equal proportions in each 
   category.
*/
   int i;
   double t, factor=alpha/beta*K, lnga1;

   if (mean) {
      lnga1=LnGamma(alpha+1);
      for (i=0; i<K-1; i++) /* cutting points, Eq. 9 */
         freqK[i]=PointGamma((i+1.0)/K, alpha, beta);
      for (i=0; i<K-1; i++) /* Eq. 10 */
         freqK[i]=IncompleteGamma(freqK[i]*beta, alpha+1, lnga1);

      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   else {
      FOR(i,K) rK[i]=PointGamma((i*2.+1)/(2.*K), alpha, beta);
      for(i=0,t=0; i<K; i++) t+=rK[i];
      FOR(i,K) rK[i]*=factor/t;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}

int MultiNomialAliasSetTable (int ncat, double prob[], double F[], int L[], double space[])
{
/* This sets up the tables F and L for the alias algorithm for generating samples from the 
   multinomial distribution MN(ncat, p) (Walker 1974; Kronmal & Peterson 1979).  
   
   F[i] has cutoff probabilities, L[i] has aliases.
   I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.

   Should perhaps check whether prob[] sums to 1.
*/
   signed char *I = (signed char *)space;
   int i,j,k, nsmall;

   for(i=0; i<ncat; i++)  L[i]=-9;
   for(i=0; i<ncat; i++)  F[i]=ncat*prob[i];
   for(i=0,nsmall=0; i<ncat; i++) {
      if(F[i]>=1)  I[i]=1;
      else       { I[i]=-1; nsmall++; }
   }
   for(i=0; nsmall>0; i++) {
      for(j=0; j<ncat; j++)  if(I[j]==-1) break;
      for(k=0; k<ncat; k++)  if(I[k]==1)  break;
      if(k==ncat)  break;

      L[j] = k;
      F[k] -= 1-F[j];
      if(F[k]<1) { I[k]=-1; nsmall++; }
      I[j]=0;  nsmall--;
   }
   return(0);
}


int MultiNomialAlias (int n, int ncat, double F[], int L[], int nobs[])
{
/* This generates multinomial samples using the F and L tables set up before, 
   using the alias algorithm (Walker 1974; Kronmal & Peterson 1979).
   
   F[i] has cutoff probabilities, L[i] has aliases.
   I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.
*/
   int i,k;
   double r;

   for(i=0; i<ncat; i++)  nobs[i]=0;
   for(i=0; i<n; i++)  {
      r = rndu()*ncat;
      k = (int)r;
      r -= k;
      if(r<=F[k]) nobs[k]++;
      else        nobs[L[k]]++;
   }
   return (0);
}     

double rndgamma1 (double s);
double rndgamma2 (double s);

double rndgamma (double s)
{
/* random standard gamma (Mean=Var=s,  with shape parameter=s, scale para=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1

   This is unsafe, and is found to return 0 when s is very small.
*/
   double r=0;

   if (s<=0)      puts ("jgl gamma..");
   else if (s<1)  r=rndgamma1 (s);
   else if (s>1)  r=rndgamma2 (s);
   else           r=-log(rndu());
   return (r);
}


double rndgamma1 (double s)
{
/* random standard gamma for s<1
   switching method
*/
   double r, x=0,small=1e-37,w;
   static double a,p,uf,ss=10,d;

   if (s!=ss) {
      a=1-s;
      p=a/(a+s*exp(-a));
      uf=p*pow(small/a,s);
      d=a*log(a);
      ss=s;
   }
   for (;;) {
      r=rndu();
      if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
      else if (r>uf)  x=a*pow(r/p,1/s), w=x;
      else            return (0);
      r=rndu ();
      if (1-r<=w && r>0)
         if (r*(w+1)>=1 || -log(r)<=w)  continue;
      break;
   }
   return (x);
}

double rndgamma2 (double s)
{
/* random standard gamma for s>1
   Best's (1978) t distribution method
*/
   double r,d,f,g,x;
   static double b,h,ss=0;
   if (s!=ss) {
      b=s-1;
      h=sqrt(3*s-0.75);
      ss=s;
   }
   for (;;) {
      r=rndu ();
      g=r-r*r;
      f=(r-0.5)*h/sqrt(g);
      x=b+f;
      if (x <= 0) continue;
      r=rndu();
      d=64*r*r*g*g*g;
      if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
   }
   return (x);
}


void ClearNode (int inode)
{
/* a source of confusion. Try not to use this routine.
*/
   nodes[inode].father=nodes[inode].ibranch=-1;
   nodes[inode].nson=0;
   nodes[inode].branch=nodes[inode].age=0;
   /* nodes[inode].label=0; */
   /* nodes[inode].branch=0; clear node structure only, not branch lengths */
   /* FOR (i, com.ns) nodes[inode].sons[i]=-1; */
}


int ReadaTreeB (FILE *ftree, int popline)
{
   char line[255];
   int nodemark[NS*2-1]={0}; /* 0: absent; 1: father only (root); 2: son */
   int i,j, state=0, YoungAncestor=0;

   if(com.clock) {
      puts("\nbranch representation of tree might not work with clock model.");
      getchar();
   }

   fscanf (ftree, "%d", &tree.nbranch);
   FOR (j, tree.nbranch) {
      FOR (i,2) {
         if (fscanf (ftree, "%d", & tree.branches[j][i]) != 1) state=-1;
         tree.branches[j][i]--;
         if(tree.branches[j][i]<0 || tree.branches[j][i]>com.ns*2-1) 
            error2("ReadaTreeB: node numbers out of range");
      }
      nodemark[tree.branches[j][1]]=2;
      if(nodemark[tree.branches[j][0]]!=2) nodemark[tree.branches[j][0]]=1;
      if (tree.branches[j][0]<com.ns)  YoungAncestor=1;

      printf ("\nBranch #%3d: %3d -> %3d",j+1,tree.branches[j][0]+1,tree.branches[j][1]+1);

   }
   if(popline) fgets(line, 254, ftree);
   for(i=0,tree.root=-1; i<tree.nbranch; i++) 
      if(nodemark[tree.branches[i][0]]!=2) tree.root=tree.branches[i][0];
   if(tree.root==-1) error2("root err");
   for(i=0; i<com.ns; i++)
      if(nodemark[i]==0) {
         matIout(F0,nodemark,1,com.ns);
         error2("branch specification of tree");
      }

   if(YoungAncestor) {
      puts("\nAncestors in the data?  Take care.");
      if(!com.cleandata) {
         puts("This kind of tree does not work with unclean data.");
         getchar();
      }
   }

/*
   com.ntime = com.clock ? (tree.nbranch+1)-com.ns+(tree.root<com.ns)
                         : tree.nbranch;
*/

   BranchToNode ();
   return (state);
}


int IsNameNumber(char line[])
{
/* returns 0 if line has species number; 1 if name; 2 if both number and name
*/
   int isname=0, j,k, ns=com.ns;
   int SeparatorFixed=(int)'_';

   if(ns<1) error2("ns=0 in IsNameNumber");
   /* both name and number? */
   k = strchr(line, SeparatorFixed) - line;
   for(j=0; j<k; j++)
      if(!isdigit(line[j])) break;
   if(j==k) 
      isname=2;
   else {
      for(j=0; line[j]; j++)  /* name or number? */
         if(!isdigit(line[j])) { isname=1; break; }  
   }
   if(isname==0 || isname==2) {
      sscanf(line,"%d",&k);
      if(k<1||k>ns) {
         printf("species number %d outside range.", k);
         exit(-1);
      }
   }
   return(isname);
}

int print1seq (FILE*fout, char *z, int ls, int encoded, int pose[])
{
/* This prints out one sequence, which may be coded.  Codon seqs are coded
   as 0,1,2,...60 if called from codeml.c or 0,1,2,...,63 otherwise.
   This may be risking.  Check when use.
   z[] contains patterns if (pose!=NULL)
   This uses com.seqtype.
*/
   int i, h,hp, gap=10;
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs)), str[4]="";
   int nb=(com.seqtype==CODONseq?3:1);


   if(!encoded) {  /* raw data, not coded */
      for(h=0; h<ls; h++) {
         hp=(pose?pose[h]:h);
         FOR(i,nb) 
            fputc(z[hp*nb+i],fout);
         if(com.seqtype==CODONseq || (h+1)%gap==0) fputc(' ',fout);
      }
   }
   else {    /* cleandata, coded */
      for(h=0; h<ls; h++) {
         hp=(pose?pose[h]:h);
         if(com.seqtype!=CODONseq) 
            fputc(pch[(int)z[hp]],fout);
         else {
#ifdef CODEML
            fprintf(fout,"%s",getcodon(str,FROM61[(int)z[hp]])); /* 0,1,...,60 */
#else
            fprintf(fout,"%s",getcodon(str,z[hp]));         /* 0,1,...,63 */
#endif
         }
         if(com.seqtype==CODONseq || (h+1)%gap==0) fputc(' ',fout);
      }
   }
   return(0);
}


long factorial(int n)
{
   long f, i;
   if (n>10) error2("n>10 in factorial");
   for (i=2,f=1; i<=(long)n; i++) f*=i;
   return (f);
}


double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   /* double accurate=1e-8, overflow=1e30; */
   double accurate=1e-10, overflow=1e60;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}

int matIout (FILE *fout, int x[], int n, int m)
{
   int i,j;
   for (i=0,FPN(fout); i<n; i++,FPN(fout)) 
      FOR(j,m) fprintf(fout,"  %4d", x[i*m+j]);
   return (0);
}


void BranchToNode (void)
{
/* tree.root need to be specified before calling this
*/
   int i, from, to;
   
   tree.nnode=tree.nbranch+1;
   for(i=0; i<tree.nnode; i++)
      { nodes[i].father=nodes[i].ibranch=-1;  nodes[i].nson=0; }
   for (i=0; i<tree.nbranch; i++) {
      from=tree.branches[i][0];
      to  =tree.branches[i][1];
      nodes[from].sons[nodes[from].nson++]=to;
      nodes[to].father=from;
      nodes[to].ibranch=i;
   }
   /*  nodes[tree.root].branch=0;  this breaks method=1 */
}


char *getcodon (char codon[], int icodon)
{
/* id : (0,63) */
   if (icodon<0||icodon>63) {
      printf("\ncodon %d\n", icodon);
      error2("getcodon.");
   }
   codon[0]=BASEs[icodon/16]; 
   codon[1]=BASEs[(icodon%16)/4];
   codon[2]=BASEs[icodon%4];
   codon[3]=0;
   return (codon);
}


double PointChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g, small=1e-6;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<small)   return(0);
   if (p>1-small) return(9999);
   if (v<=0)      return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3: 
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0)
      error2 ("\nIncompleteGamma");
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}

double PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.
*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) z=999;
   else {
      y = sqrt (log(1/(p1*p1)));   
      z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   }
   return (p<0.5 ? -z : z);
}

#endif /* _MCMCCOAL_C_ */
