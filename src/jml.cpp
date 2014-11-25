/**********************************************\
 |                                              |
 |                     JML                      |
 |                                              |
 |               version 1.3.0                  |
 |                                              |
 |----------------------------------------------|
 |                                              |
 |                                              |
 |        released November 24th, 2014          |
 |                                              |
 |   (c) Copyright 2011-2014 by Simon Joly      |
 |                                              |
 \**********************************************/

/*
	JML, version 1.3.0
    copyright (c) Simon Joly, 2011-2014
		
    JML is a free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2.

    JML is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include "jml.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <cctype>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include "treeinfo.h"                   // define classes for parsing newick trees
#include "nexusdata.h"

#define LARGENUMBER 10000
#define isNL_NEXUSwhiteSpace(c)  ( strchr(" \t\v\f", (c)) || (((c) <= 6) && ((c) >= 0)))

using namespace std;

void outputusage(void);
void readcontrolfile(void);
void convertseqgencommands(void);
string To_Uppercase(string input_string);
void OutputTimeUsed(void);
template<typename T>
extern string to_string(const T & Value);
extern double string_to_double(string a_string);
void CalMinDist(void);
int string_to_int(string a_string)
{
	int number;
	istringstream iss;
	iss.str(a_string);
	iss >> number;
	return(number);
}

nexusdata nexus_data, obs_data;
string input_file = "species.trees";
string control_file = "jml.ctl";
string data_file = "";
ofstream outfile, batchfile, histo, mindistmat;
vector<string> speciesnames;
vector<int> seqinspecies;
int numberofspecies;
string seqgencommand="";
int seqgen_argc;
char **seqgen_argv;
double locusrate = 1;
double heredityscalar = 1;
double significancelevel = 0.1;
int burnin=-1;
time_t start_time;
int seconds;
vector < vector <double> > mindist;  //Array containint the posterior distributions...
string tax1,tax2;
int pos1,pos2;
double tempdist,smallestdist;
string temp;
int seed=-1;  //seed
int jump=-1;  //parameter to store the thinning value


int main(int argc, char **argv)
{
	char *cur;
	char ch;
	int i,w,c;
	vector<Tree> all_trees;
	all_trees.reserve(200000);

	cout << endl;
	cout << " ************************************************************************" << endl;
	cout << " *                                                                      *" << endl;
	cout << " *  JML, version " << VERSION;
	cout <<                     "                                                    *" << endl;
	cout << " *  Copyright 2011-2014 by Simon Joly                                   *" << endl;
	cout << " *                                                                      *" << endl;
	cout << " *  ------------------------------------------------------------------  *" << endl;
	cout << " *                                                                      *" << endl;
	cout << " *  This is a free software and it comes with absolutely no warranty.   *" << endl;
	cout << " *  You can redistribute the software and/or modify it under the terms  *" << endl;
	cout << " *  of the GNU General Public License version 3. See the GNU General    *" << endl;
	cout << " *  Public License for more details (file gpl.txt that comes with the   *" << endl;
	cout << " *  source files).                                                      *" << endl;
	cout << " *                                                                      *" << endl;
	cout << " ************************************************************************" << endl;
	cout << endl;

    if(argc<2) {
        outputusage();
        exit(0);
    }
	
	for(w=1;w<argc;w++)
		{
		cur=argv[w];
		ch=cur[0];
	
		if (ch == '-')
			{
			ch=toupper(cur[1]);
			switch(ch)
				{
				case 'T':
					w++;
					if(w<argc && argv[w][0] != '-') 
						{
						input_file = argv[w];
						}
					else {w--;}	
					continue;
                case 'H':
                    outputusage();
                    exit(0);
				case 'C':
					w++;
					if(w<argc && argv[w][0] != '-')    
						{
						control_file = argv[w];
						}
					else {w--;}	
					continue;
				case 'D':
					w++;
					if(w<argc && argv[w][0] != '-')    
						{
						data_file = argv[w];
						}
					else {w--;}	
					continue;
                case 'S':
                    w++;
                    if(w<argc && argv[w][0] != '-')    
						{
                            data_file = argv[w];
						}
                    else {w--;}	
                    continue;
				default:
						continue;
				}
			}
		}


	readcontrolfile();
	convertseqgencommands();
	if (data_file != "") CalMinDist();
	Tree a_tree;
	a_tree.enterlocusrate(locusrate);
	a_tree.enterheredityscalar(heredityscalar);
	vector<string> taxnumbers;
	vector<string> taxlabels;

	ifstream infile;
	for (;;) {
		infile.open(input_file.c_str());
		if (!infile) {
				cout << " Error in opening species tree file: " << input_file << endl;
				cout << " Enter the correct name of the species tree file: ";
				cin >> input_file;
				continue;
				}
		break;
		}        
	cout << endl << " Reading tree file \"" << input_file << "\"... ";
    cout.flush();
	infile >> temp;
	while(1) {
		while (To_Uppercase(temp) != "BEGIN") infile >> temp;
		infile >> temp;
		if (To_Uppercase(temp) == "TREES;") {
			infile >> temp;
			if (To_Uppercase(temp) == "TRANSLATE") {  // Translare block present
				a_tree.translateblock();
				infile >> temp;
				while (temp != ";") {
					taxnumbers.push_back(temp);
					infile >> temp;
					int c = (int)temp.length()-1;
					while(strchr(" ,",temp[c])) {
						temp.erase(c);
						c--;
						}
					taxlabels.push_back(temp);
					infile >> temp;
					}
				a_tree.gettaxnumbers(taxnumbers);
				a_tree.gettaxlabels(taxlabels);
				}
			while(1) {
				getline(infile,temp);
				if (To_Uppercase(temp) == "END;") break;
				else if (temp == "") continue;
				else {
					a_tree.parseTree(temp);
					all_trees.push_back(a_tree);
					}
				}
			break;
			}
		else continue;
		}
    infile.close();
    cout << "done." << endl;

    if (burnin == -1){
        burnin = ceil(0.1*all_trees.size());
    }
	if (jump < 1) jump = 1;

	/* initialize random seed: */
    if (seed==-1) {
        seed = time(NULL);
        srand ( seed );
    }
	else srand(seed);
	unsigned long  randomseed;
	
	int trees_post_burnin = (all_trees.size() - burnin);
	int sampled_trees = (int)(trees_post_burnin/jump);

	cout << endl;
	cout << " *** Species trees file information ***" << endl;
	cout << " Number of trees read: " << all_trees.size() << endl;
	cout << " Number of trees in burnin: " << burnin << endl;
	cout << " Number of trees post-burnin: " << trees_post_burnin << endl;
    //cout << " Thinning: " << jump << endl;
	cout << " Number of trees sampled for simulations: " << sampled_trees << endl;
    cout << " Seed used: " << seed << endl;
	cout << " Significance level: " << significancelevel << endl << endl;
	cout << " *** Locus-specific information ***" << endl;
	cout << " Locus relative mutation rate: " << locusrate << endl;
	cout << " Locus heredity scalar: " << heredityscalar << endl << endl;

	int z;
	int j,k,l,m;
	
	start_time=time(NULL);
	
	//Array containint the posterior distributions...
	vector< vector < vector <double> > > postdist;
	postdist.resize(numberofspecies);
	for(k=0;k<numberofspecies;k++) {
		postdist[k].resize(numberofspecies);
		}
	for(i=0;i<sampled_trees;i++) {
		printf("\r Performing simulation %i",(i+1));
		fflush(stdout);
		randomseed = rand();
		outfile.open("Rep.dat",ios_base::trunc);
		if(!outfile) {cout << "error in opening outfile" << endl; exit(0);}
		outfile << "Rep.phy" << endl;
		outfile << randomseed << endl;
		outfile << numberofspecies;
		for (w=0;w<numberofspecies;w++)
			outfile << " " << speciesnames[w];
		outfile << endl;
		for (w=0;w<numberofspecies;w++)
			outfile << seqinspecies[w] << " ";
		outfile << endl;
		if (i==0) outfile << all_trees[i].PrintMCcoalTree();
		else outfile << all_trees[(i*jump)].PrintMCcoalTree();
		outfile.flush();
		outfile.close();
		
		// Run MCcoal simul
		mcmccoal((char *)"Rep.dat");

		// Run seq-gen
		seqgen(seqgen_argc, seqgen_argv,(char *) "Rep.phy",randomseed);
		//cout << "  Out seq-gen" << endl;

		//Analyse simulated data
		infile.open("Rep.phy");
		if (!infile) {cerr << "Error in opening file Rep.phy" << endl;exit(0);} 
		
		/*** Read Phylip file ***/
		nexus_data.InitializeMember();
		infile >> temp;
		nexus_data.GetNTaxa(atoi(temp.c_str()));
		infile >> temp;
		nexus_data.GetNChars(atoi(temp.c_str()));
		nexus_data.InitDataMatrix();
		nexus_data.InitMatrix();
		for (m=0;m<nexus_data.ReturnNTax();m++)
				{
				infile >> temp;
				nexus_data.ReadTaxa(temp);
				while ((c=infile.get()) != '\n')
						{
						if (isNL_NEXUSwhiteSpace(c)) continue;
						temp = c;
						nexus_data.AddCharacter(m, temp);
						}
				}
		infile.close();
		//cout << "  Phylip file read" << endl;
	
		// Compute distances between pairs of sequences
		nexus_data.computedistances();
		
		for(k=0;k<numberofspecies;k++) {
			for(l=0;l<k;l++) {
				smallestdist=LARGENUMBER;
				for (m=0;m<seqinspecies[k];m++) {
					for (j=0;j<seqinspecies[l];j++) {
// [SJ] changes in version 1.03: there is now always a number after the species names,
//        even when there is a single sequence for that species.
//						if (seqinspecies[k] == 1) tax1 = speciesnames[k]; // Only one seq in species k
//						else tax1 = speciesnames[k] + to_string(m+1);
						tax1 = speciesnames[k] + to_string(m+1);
//						if (seqinspecies[l] == 1) tax1 = speciesnames[l]; // Only one seq in species l
//						else tax2 = speciesnames[l] + to_string(j+1);
						tax2 = speciesnames[l] + to_string(j+1);
						pos1 = nexus_data.GetPositionofTaxa(tax1); pos2 = nexus_data.GetPositionofTaxa(tax2);
						tempdist = nexus_data.ReturnDist(pos1,pos2);
						if (tempdist < smallestdist) smallestdist = tempdist;
						}
					}
				postdist[k][l].push_back(smallestdist);
				}
			}
		//cout << "one simul done" << endl;			
		//getchar();
		}
	
	cout << endl << endl << " Writing Posterior Predictive Distributions" << endl;

	//Sort vectors
	for(k=0;k<numberofspecies;k++) {
		for(l=0;l<k;l++) {
			sort(postdist[k][l].begin(), postdist[k][l].end());
			}
		}

	// Print a file that contains posterior predictive values for all species comparisons
	histo.open("Distributions.txt");
	if (!histo) {
			cout << "Error in openning batchfile" << endl;
			exit(0);
			}
	// column names
	for(k=0;k<numberofspecies;k++) {
		for(l=0;l<k;l++) {
			histo << (speciesnames[k]+"-"+speciesnames[l]) << '\t' ;
			}
		}
	histo << endl;
	for(z=0;z<sampled_trees;z++)
		{
		for(k=0;k<numberofspecies;k++) {
			for(l=0;l<k;l++) {
				histo << postdist[k][l][z] << '\t';
				}
			}
		histo << endl;
		}
	histo.close();

	if(data_file != "") {
		//Calculate probabilities
		cout << " Calculating probabilities..." << endl;
		ofstream probfile;
		probfile.open("Probabilities.txt");
		if (!probfile) {
				cout << "Error in openning Probabilities file" << endl;
				exit(0);
				}
		vector < vector <double> > Probs;
		Probs.resize(numberofspecies);
		for(k=0;k<numberofspecies;k++) {
			Probs[k].resize(numberofspecies);
			}
		for(k=0;k<numberofspecies;k++) {
			for(l=0;l<numberofspecies;l++) {
				Probs[k][l] = 1;
				}
			}
		for(k=0;k<numberofspecies;k++) {
			for(l=0;l<k;l++) {
				for(m=0;m<sampled_trees;m++) {
					if (postdist[k][l][m]>mindist[k][l]) {
						Probs[k][l] = ((double)(m+1)/sampled_trees);
						break;
						}
					}
				}
			}
		//Print probabilities
		probfile << "Comparison\tminDist\tProbability" << endl;
		for(k=0;k<numberofspecies;k++) {
			for(l=0;l<k;l++) {
				probfile << (speciesnames[k]+"-"+speciesnames[l]) << '\t' << mindist[k][l] << '\t' << Probs[k][l] << endl;
				}
			}
		probfile.close();

        //Calculate probabilities for all distances with p < significancelevel
        ofstream resultsfile;
        resultsfile.open("Results.txt");
        if (!resultsfile) {
                cout << "Error in openning Results file" << endl;
                exit(0);
                }
        resultsfile << "Sp Comparison\tseq1\tseq2\tDistance\tProbability" << endl;
        int mm;
        for(k=0;k<numberofspecies;k++) {
            for(l=0;l<k;l++) {
                smallestdist=LARGENUMBER;
                for (m=0;m<seqinspecies[k];m++) {
                    for (j=0;j<seqinspecies[l];j++) {
                        tax1 = speciesnames[k] + to_string(m+1);
                        tax2 = speciesnames[l] + to_string(j+1);
                        pos1 = obs_data.GetPositionofTaxa(tax1); pos2 = obs_data.GetPositionofTaxa(tax2);
                        tempdist = obs_data.ReturnDist(pos1,pos2);
                        if (tempdist < postdist[k][l][int(significancelevel*sampled_trees)]) {
                            //Calculate exact probability...
                            for(mm=0;mm<sampled_trees;mm++) {
                                if (postdist[k][l][mm]>tempdist) {
                                    resultsfile << (speciesnames[k]+"-"+speciesnames[l]) << '\t' << tax1 << '\t' << tax2 << '\t';
                                    resultsfile << tempdist << '\t' << ((double)(mm+1)/sampled_trees) << endl;
                                    break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        resultsfile.close();
        }
		
	cout << endl << " One example tree (check to make sure it looks ok):" << endl;
	cout << endl << a_tree.PrintMCcoalTree() << endl << endl;
	cout << " Simulations done." << endl;
	OutputTimeUsed();
	cout << " Posterior predictive distributions can be found in file \"Distributions.txt\"" << endl;
	if(data_file != "") {
		cout << " Probabilities of the minimum distances can be found in the file \"Probabilities.txt\"" << endl;
		//cout << " The matrix of minimum distances can be found in the file \"MinimumDistances.txt\"" << endl;
		cout << " Significant results (p < " << significancelevel << ") can be found in the file \"Results.txt\"" << endl;
		}
    cout << endl;
	return 0;
}


void CalMinDist(void)
{
    cout << " Calculating minimum distances from file... ";
    cout.flush();
	static int k,l,m,c,j; 
	ifstream datafile;
	mindist.resize(numberofspecies);
	for(k=0;k<numberofspecies;k++) {
		mindist[k].resize(numberofspecies);
		}
	//Analyse simulated data
	for (;;) {
		datafile.open(data_file.c_str());
		if (!datafile) {
				cout << " Error in opening data file: " << data_file << endl;
				cout << " Enter the correct name of the sequence data file: ";
				cin >> data_file;
				continue;
				}
		break;
		}        

	/*** Read Phylip file ***/
	//nexusdata obs_data;
	obs_data.InitializeMember();
	datafile >> temp;
	obs_data.GetNTaxa(atoi(temp.c_str()));
	datafile >> temp;
	obs_data.GetNChars(atoi(temp.c_str()));
	obs_data.InitDataMatrix();
	obs_data.InitMatrix();
	for (m=0;m<obs_data.ReturnNTax();m++)
			{
			datafile >> temp;
			obs_data.ReadTaxa(temp);
			while ((c=datafile.get()) != '\n')
					{
					if (isNL_NEXUSwhiteSpace(c)) continue;
					temp = c;
					obs_data.AddCharacter(m, temp);
					}
			}
	datafile.close();
	
	// Calculate minimum distances between pairs of sequences
	obs_data.computedistances();

	for(k=0;k<numberofspecies;k++) {
		for(l=0;l<k;l++) {
			smallestdist=LARGENUMBER;
			for (m=0;m<seqinspecies[k];m++) {
				for (j=0;j<seqinspecies[l];j++) {
					tax1 = speciesnames[k] + to_string(m+1);
					tax2 = speciesnames[l] + to_string(j+1);
					pos1 = obs_data.GetPositionofTaxa(tax1); pos2 = obs_data.GetPositionofTaxa(tax2);
					tempdist = obs_data.ReturnDist(pos1,pos2);
					if (tempdist < smallestdist) smallestdist = tempdist;
					}
				}
			mindist[k][l] = smallestdist;
			}
		}
    cout << "done." << endl;
		
	//Output minimum distances to file
/*	mindistmat.open("MinimumDistances.txt");
	if (!mindistmat) {
			cout << "Error in openning MinimumDistances file" << endl;
			exit(0);
			}
	for(k=0;k<numberofspecies;k++) mindistmat  << '\t' << speciesnames[k];
	mindistmat << endl;
	for(k=0;k<numberofspecies;k++) {
		mindistmat << speciesnames[k];
		for(l=0;l<k;l++) {
			mindistmat << '\t' << mindist[k][l];
			}
		mindistmat << endl;
		}
	mindistmat.close();
*/
}

void outputusage(void)
{
    cout << endl;
    cout << " Usage for JML:" << endl << endl;
    cout << "   flag  description              Default name   Required?" << endl;
    cout << "   ----  ----------------------   -------------  ---------" << endl;
    cout << "    -c   <control_file_name>      jml.ctl        yes" << endl;
    cout << "    -t   <input_tree_file_name>   species.trees  yes" << endl;
    cout << "   [-d]  <sequence_file_name>     [none]         no" << endl;
    cout << endl << endl;
    cout << " Example:" << endl << endl;
    cout << "   [UNIX]  ./jml  -c jml.ctl  -t species.trees  -d mysequences.phy" << endl;
    cout << "   [PC]    jml  -c jml.ctl  -t species.trees  -d mysequences.phy" << endl;
    cout << endl << endl;
}


void OutputTimeUsed(void)
{
	time_t now=time(NULL);
	seconds = now - start_time;
	cout << " Time used (Hr:Min:Sec): ";
	int prev_width=cout.width(2);
	char prev=cout.fill('0');
	cout << right << floor(seconds/3600);
	cout.width(prev_width);cout.fill(prev);
	cout << ":";
	cout.width(2);cout.fill('0');
	cout << floor( (seconds - (floor(seconds/3600)*3600))/60 );
	cout.fill(prev);cout.width(prev_width);
	cout << ":";
	cout.width(2);cout.fill('0');
	cout << (seconds - ( (floor(seconds/3600)*3600) + (floor( (seconds - (floor(seconds/3600)*3600))/60)*60) ) ) << endl;
	cout.fill(prev);cout.width(prev_width);
}



/*

Function: readcontrolfile

Detail: Open the control file (jml.ctl by default) and reads its content

*/

void readcontrolfile(void)
{
    cout << " Reading control file... ";
    cout.flush();
	string lineofdata;
	string tempstring;
	static int cursor;
	ifstream controlfile;
	for (;;) {
		controlfile.open(control_file.c_str());
		if (!controlfile) {
				cout << "Error in opening control file: " << control_file << endl;
				cout << "Enter the correct control file name: ";
				cin >> control_file;
				continue;
				}
		break;
		}        
	
	while(!controlfile.eof()) {
		getline(controlfile,lineofdata);
//cout << lineofdata << endl;
		cursor=0;
		tempstring="";
		if ( (cursor) == (int)lineofdata.length()) break;       //Skip white lines...
		while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
		if ( (cursor) == (int)lineofdata.length()) break;       //Skip white lines...
		while (!strchr("= \t\v\r\f\n",lineofdata[cursor])) {
			tempstring += lineofdata[cursor];
			cursor++;
			}
		while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
		if (lineofdata[cursor] != '=') cerr << "Command " << tempstring << " is not followed by an \"=\" sign" << endl;
		cursor++;
		while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
		if (To_Uppercase(tempstring) == "SPECIES") {
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
					}
				speciesnames.push_back(tempstring);
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor]) && (cursor != (int)lineofdata.length())) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
				}
			}
		else if (To_Uppercase(tempstring) == "SEQPERSPECIES") {
//cout << endl << "Sequences per species" << endl;
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
					}
//cout << tempstring << " ";
				seqinspecies.push_back(string_to_int(tempstring));
//cout << string_to_int(tempstring) << " ";
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor]) && (cursor != (int)lineofdata.length())) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
				}			
			}
		else if (To_Uppercase(tempstring) == "LOCUSRATE") {
//cout << endl << "Locus rate" << endl;
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
					}
//cout << tempstring << " ";
				locusrate = string_to_double(tempstring);
//cout << string_to_int(tempstring) << " ";
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
				}			
			}
		else if (To_Uppercase(tempstring) == "HEREDITYSCALAR") {
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
					}
				heredityscalar = string_to_double(tempstring);
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
				}
			}
		else if (To_Uppercase(tempstring) == "SEQGENCOMMAND") {
			while (1) {
				seqgencommand="";
				while (cursor != (int)lineofdata.length() ) {
					seqgencommand += lineofdata[cursor];
					cursor++;
					}
				break;
				}
			string::iterator it=seqgencommand.end();
			int b = (int)seqgencommand.length()-1;
			while(strchr(" \t\v\r\f\n",seqgencommand[b])) {
				seqgencommand.erase(it);
				it--;b--;
				}
			}
		else if (To_Uppercase(tempstring) == "SIGNIFICANCELEVEL") {
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
                }
				significancelevel = string_to_double(tempstring);
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
            }
        }
		else if (To_Uppercase(tempstring) == "BURNIN") {
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
                }
				burnin = string_to_int(tempstring);
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
            }
        }
		else if (To_Uppercase(tempstring) == "THINNING") {
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
                }
				jump = string_to_int(tempstring);
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
            }
        }
		else if (To_Uppercase(tempstring) == "SEED") {
			while (1) {
				tempstring="";
				while ((!strchr(" \t\v\r\f\n",lineofdata[cursor])) && (cursor != (int)lineofdata.length()) ) {
					tempstring += lineofdata[cursor];
					cursor++;
                }
				seed = string_to_int(tempstring);
				if ( (cursor) == (int)lineofdata.length()) break;
				while (strchr(" \t\v\r\f\n",lineofdata[cursor])) cursor++;  //strip white spaces
				if ( (cursor) == (int)lineofdata.length()) break;
            }
        }
    }
    cout << "done." << endl;
	
	numberofspecies = (int)speciesnames.size();
	// Few test is data was entered correctly
	if (speciesnames.size() ==0) {
		cerr << " No species information entered... Exiting program!" << endl;
		exit(0);
		}
	if (seqinspecies.size() ==0) {
		cerr << " No sequence number information entered... Exiting program!" << endl;
		exit(0);
		}
	if (seqinspecies.size() != speciesnames.size()) {
		cerr << " The number of species (" << speciesnames.size() << ") is not equal to the number of sequences per species (";
		cerr << seqinspecies.size() << ") provided... Exiting program!" << endl;
		exit(0);
		}
	if (seqgencommand=="") {
		cerr << " No seq-gen commands entered... Exiting program!" << endl;
		exit(0);
		}
		
}

/* Convert the seq-gen command line in an array that will then be passed to seq-gen for simulations */

void convertseqgencommands(void)
{
	static int thiscursor=0;
	static int h,g;
	string temporary_string="";
	vector<string> arguments;
	for(;;) {
		if (thiscursor == (int)seqgencommand.length()) break;
		if (strchr(" \t\v\f\n\r",seqgencommand[thiscursor])) { // character is a white space
			thiscursor++;
			continue;
			}
		while ((!strchr(" \t\v\f\n\r",seqgencommand[thiscursor])) && (thiscursor != (int)seqgencommand.length()) ) {
			temporary_string += seqgencommand[thiscursor];
			thiscursor++;
			}
		arguments.push_back(temporary_string);
		temporary_string = "";
		//if (thiscursor == (int)seqgencommand.length()) break;
		}
	seqgen_argv = new char* [arguments.size()+2];
	seqgen_argv[0] = (char *) "-q"; // quiet
	seqgen_argv[1] = (char *) "out.trees"; // tree file
	for (h=2;h<((int)arguments.size()+2);h++) {
		seqgen_argv[h] = new char[arguments[(h-2)].size()];
		for (g=0;g<(int)arguments[(h-2)].size();g++) {
			seqgen_argv[h][g] = arguments[(h-2)][g];
			}
		seqgen_argv[h][arguments[(h-2)].size()] = '\0';
		}
	seqgen_argc = (int)(arguments.size()+2);
	/*
	cout << endl << endl;
	for (h=0;h<seqgen_argc;h++) {
		cout << seqgen_argv[h] << " ";
		}
	*/
}


string To_Uppercase(string input_string)
{
    int i;
    char character;
    string output_string, temp_string;
    temp_string = input_string;
    int size_of_string = temp_string.size();

    for (i=0; i < size_of_string; i++)
    {
        character = temp_string[i];
        output_string += toupper(character);
    }
    return output_string;
}
