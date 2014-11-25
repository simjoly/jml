/* MCMCcoal.h
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


#ifndef _MCMCCOAL_H_
#define _MCMCCOAL_H_

#ifdef _cplusplus_ /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

void mcmccoal(char *simul_file_name);

#ifdef _cplusplus_ /* If this is a C++ compiler, end C linkage */
}
#endif

#endif /* _MCMCCOAL_H_ */
