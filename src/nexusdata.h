/***************************\
*                           *
*      nexusdata class      *
*                           *
\***************************/

/*

		Class Nexusdata
    (c) Copyright 2005-2011 by Simon Joly

    This file is part of JML.

    JML is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Pofad is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <string>
#include <iostream>
#include <vector>
#include "GenFunctions.h"

using namespace std;

#ifndef _NEXUSDATA_CLASS_DEFINED_

class nexusdata {
    private:
        int NTaxa;            // number of taxa in the dataset
        int NChars;           // The number of characters
        vector<string> Taxa;  // A vector containning the taxa's name
        int TaxaRead;         // A flag to show if taxa labels are read: 0=not read; 1=read.
        int diagonal;         // 0 = no diagonal, 1 = diagonal
		int Labels;           // Flag to indicate is labels are in the distance matrix: 0=NoLabels, 1=Labels
        int interleave;       // Flag to indicate if matrix is interleaved: 0=no; 1=yes
        string datatype;      // Indicate the dataype
        string missingchar;   // Character indicating the symbol for missing characters
        string gapchar;       // Character for gaps
        int triangle;         // 0=Lower, 1=Upper, 2=both
        double **dist_matrix; // The distance matrix
        string *char_matrix;  // The character matrix

    public:
        nexusdata();
		~nexusdata();
        void InitializeMember();
        void GetNTaxa(int nb_taxa);
        int ReturnNTax();
        void GetNChars(int nb_chars);
        int ReturnNChars();
        void ReadTaxa(string a_taxa);
		string ReturnTaxa(int number);
        int NumberTaxaLabels();
        void isTaxa();
        void IsDiagonal(int diag);
        int ReturnDiagonal();
        void isLabels(int number);
        int ReturnIsLabels();
        void IsInterleave(int a_number);
        int ReturnInterleave();
        void MissingChar(string a_character);
        string ReturnMissingChar();
        void GetGapChar(string a_character);
        string ReturnGapChar();
        void InitMatrix();
        void IsTriangle(int a_number);
        int ReturnTriangle(void);
        void EnterMatrix(double distance, int i, int j);
        double ReturnDist(int i, int j);

        int GetPositionofTaxa(string a_taxa);
        int ReturnisTaxa();
        void GetDatatype(string the_datatype);
        string ReturnDatatype();
        void InitDataMatrix();
        void AddCharacter(int i, string characters);
        int ReturnCharacterForTaxa(int i);
        char ReturnChar(int i, int j);                    //Returns character at position j+1 for Taxa[i]
        string ReturnSequence(int i);                     //Returns the sequence of Taxa[i]
        void computedistances(void);
};


inline nexusdata::nexusdata() {
    TaxaRead = 0;
    diagonal = 1;
    triangle = 0;
    Labels = 1;
    interleave = 0;
    }
inline nexusdata::~nexusdata(){
		int f;
		for (f=0; f < NTaxa; f++)
			delete dist_matrix[f];
    delete [] dist_matrix;
    delete [] char_matrix;
    }
inline void nexusdata::InitializeMember(){
		int f;
		for (f=0; f < NTaxa; f++)
			delete dist_matrix[f];
    delete [] dist_matrix;
    delete [] char_matrix;
    NTaxa = 0;
    NChars = 0;
    Taxa.clear();
    TaxaRead = 0;
    diagonal = 1;
    Labels = 1;
    interleave = 0;
    datatype.clear();
    missingchar.clear();
    gapchar.clear();
    triangle = 0;
    }
inline void nexusdata::GetNTaxa(int nb_taxa) {
    this->NTaxa = nb_taxa;
    }
inline int nexusdata::ReturnNTax() {
    return NTaxa;
    }
inline void nexusdata::GetNChars(int nb_chars) {
    this->NChars = nb_chars;
    }
inline int nexusdata::ReturnNChars() {
    return NChars;
    }
inline void nexusdata::ReadTaxa(string a_taxa) {
    Taxa.push_back(a_taxa);
    }
inline string nexusdata::ReturnTaxa(int number) {
    return Taxa[number];
    }
inline int nexusdata::NumberTaxaLabels() {
    return Taxa.size();
    }
inline void nexusdata::isTaxa() {
    this->TaxaRead = 1;
    }
inline void nexusdata::IsDiagonal(int diag) {
    this->diagonal = diag;
    }
inline int nexusdata::ReturnDiagonal() {
    return diagonal;
    }
inline void nexusdata::isLabels(int number) {
    this->Labels = number;
    }
inline int nexusdata::ReturnIsLabels() {
    return Labels;
    }
inline void nexusdata::IsInterleave(int a_number) {
    this->interleave = a_number;
    }
inline int nexusdata::ReturnInterleave() {
    return interleave;
    }
inline void nexusdata::MissingChar(string a_character) {
    this->missingchar = a_character;
    }
inline string nexusdata::ReturnMissingChar() {
    return missingchar;
    }
inline void nexusdata::GetGapChar(string a_character) {
    this->gapchar = a_character;
    }
inline string nexusdata::ReturnGapChar() {
    return gapchar;
    }
inline void nexusdata::InitMatrix() {
    int x,y;
    dist_matrix = new double* [NTaxa];   //Allocate place for matrix
    for (y=0; y < NTaxa; y++)
        dist_matrix[y] = new double[NTaxa];
    if (!dist_matrix)
        {
        cerr << "Can't allocate space for creating the distance matrix of haplotype" << endl;
        }
    for (x=0; x < NTaxa; x++) //Initialize the matrix with 0s
        {
        for (y=0; y < NTaxa; y++)
            {
            dist_matrix[x][y] = 0;
            }
        }
    }
inline void nexusdata::IsTriangle(int a_number) {
    this->triangle = a_number;
    }
inline int nexusdata::ReturnTriangle(void){
    return triangle;
    }
inline void nexusdata::EnterMatrix(double distance, int i, int j) {
    dist_matrix[i][j] = distance;
    }
inline double nexusdata::ReturnDist(int i, int j) {

    #ifdef DEBUG
    cout << "nex_data:" << i << "," << j << ";";
    cout << dist_matrix[i][j];
    #endif

    return dist_matrix[i][j];
    }

/**** New Functions ****/

inline int nexusdata::GetPositionofTaxa(string a_taxa) {
    unsigned int i;
    int position;
    vector<string> a_vector;   //TODO: essayer a_vector = Taxa
    for (i=0; i < Taxa.size() ; i++)
        {
        a_vector.push_back(Taxa[i]);
        }
    vector<string>::iterator an_iterator;
    an_iterator = find(a_vector.begin(), a_vector.end(), a_taxa);
    if (an_iterator == a_vector.end())
        {
        cout << "The taxa " << a_taxa << " was not found in the matrix!" << endl;
        exit(1);
        }
    else
        {
        a_vector.erase(an_iterator, a_vector.end());
        position = a_vector.size();
        } 
    return position;
    }

inline int nexusdata::ReturnisTaxa() {
    return TaxaRead;
    }

inline void nexusdata::GetDatatype(string the_datatype) {
    this->datatype = the_datatype;
    }
inline string nexusdata::ReturnDatatype() {
    return datatype;
    }
inline void nexusdata::InitDataMatrix() {
//    delete [] char_matrix;
    char_matrix = new string[NTaxa];
    if (!char_matrix) {
        cerr << "Can't allocate space for creating the character matrix" << endl;
        }
    int i;
    /* Initialise the matrix */
    for (i=0; i<NTaxa; i++) {
        char_matrix[i] = "";
        }
    }
inline void nexusdata::AddCharacter(int i, string characters) {
    char_matrix[i] += characters;
    }
inline int nexusdata::ReturnCharacterForTaxa(int i) {
    return char_matrix[i].size();
    }
inline char nexusdata::ReturnChar(int i, int j) {
    return char_matrix[i][j];
    }
inline string nexusdata::ReturnSequence(int i) {
    return char_matrix[i];
    }
inline void nexusdata::computedistances(void) {
    int i,j,k;
    double final_distance,total_distance;
    int numberchar;

    for (i=0;i<NTaxa;i++) {
        for (j=0;j<=i;j++) {
            total_distance = 0;
            final_distance = 0;
            numberchar = 0;
            //#pragma omp parallel for schedule(dynamic)
            for(k=0;k<NChars;k++) {
				if ( (char_matrix[i][k] == '-') || (char_matrix[j][k] == '-') ) continue; //ignore gaps
                if (char_matrix[i][k] != char_matrix[j][k]) total_distance += 1;
                numberchar++;
                }
            final_distance = (total_distance / NChars);
            dist_matrix[i][j] = final_distance;
            dist_matrix[j][i] = dist_matrix[i][j];
            //cout << i << " " << j << " -> Dist: " << final_distance << endl;
            }
        }
    }


#define _NEXUSDATA_CLASS_DEFINED_
#endif /* _NEXUSDATA_CLASS_DEFINED_ */
