/***************************\
*                           *
*      treeinfo class       *
*                           *
\***************************/

/*

		Class treeinfo
    (c) Copyright 2011 by Simon Joly

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
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

template<typename T>
string to_string( const T & Value )
{
    // utiliser un flux de sortie pour créer la chaîne
    ostringstream oss;
    // Écrire la valeur dans le flux
    oss.precision(12);
    oss << Value;
    // renvoyer une string
    return oss.str();
}

double string_to_double(string a_string)
{
	double number;
	// créer un flux à partir de la chaîne à convertir
	istringstream iss;
	iss.str(a_string);
	// convertir en un double
	iss >> number;
	//cout << setprecision(16) << number << endl;
	return(number);
}

#ifndef _Tree_CLASS_DEFINED_

// reference class to parse the string containing the trees
class Reference {
	public:
		int count;
		int cursor;
		Reference()
			{
			count=0;
			cursor=0;
			}
};

class Node {
	public:
		Node* leftdesc;
		Node* rightdesc;
		Node* ancestor;
		double brlength;
		double popsize;
		string species;
		bool tip;
		int number_of_visits; //Used when printing trees
		bool isroot;
		Node() {
			leftdesc = NULL;
			rightdesc = NULL;
			ancestor = NULL;
			isroot = false;
			brlength = -999;
			popsize = -999;
			tip=false;
			number_of_visits=0;
			}
};


class Tree {
	private:
		Node* root;
		string newicktree;  // used to store the tree for printing
		int cursor;
		string tree_string; //used when parsing the tree
		bool translate;
		vector<string> taxnumbers;
		vector<string> taxlabels;
		double heredityscalar;
		double locusrate;
	
  public:
    Tree() {
			root=NULL;
			cursor=0;
			translate=false;
			heredityscalar=1.0;
			locusrate=1.0;
			}
		~Tree(){
			//delete root;
			}
		string PrintTree(void);
		string PrintMCcoalTree(void);
		void Print_Tree(Node*);
		void Print_MCcoal_Tree(Node*);
		void parseTree(string);
		void nextnode(Node*& actual_node);
		double getnodeage(Node* p);
		void translateblock(void);
		void gettaxnumbers(vector<string> a_vector);
		void gettaxlabels(vector<string> a_vector);
		void enterheredityscalar(double a_number);
		void enterlocusrate(double a_number);
		
};

string Tree::PrintTree(void)
	{
	newicktree = "";
	Print_Tree(root);
	return(newicktree);
	}

void Tree::enterheredityscalar(double a_number)
	{
	this->heredityscalar = a_number;
	}

void Tree::enterlocusrate(double a_number)
	{
	this->locusrate = a_number;
	}
	
void Tree::translateblock(void)
	{
	this->translate = true;
	}

void Tree::gettaxnumbers(vector<string> a_vector)
	{
	this->taxnumbers = a_vector;
	}

void Tree::gettaxlabels(vector<string> a_vector)
	{
	this->taxlabels = a_vector;
	}
	
string Tree::PrintMCcoalTree(void)
	{
	newicktree = "";
	Print_MCcoal_Tree(root);
	return(newicktree);
	}
	
double Tree::getnodeage(Node* p)
{
	if (p->leftdesc->tip == true) {
		return p->leftdesc->brlength;   // no descendant nodes
		}
	else return p->leftdesc->brlength + getnodeage(p->leftdesc);
}

void Tree::Print_MCcoal_Tree(Node* p)
	{
	if (p->tip) {
		if(this->translate) newicktree+=taxlabels[(int)string_to_double(p->species)-1];  // There was a translate block
		else newicktree+=p->species;
		if (p->popsize != -999) {
			newicktree+=(" #"+to_string((p->popsize*this->heredityscalar*this->locusrate)));
			}
		Print_MCcoal_Tree(p->ancestor);
		}
	else {
		if (p->number_of_visits == 0) {
			p->number_of_visits++;
			newicktree+="(";
			Print_MCcoal_Tree(p->leftdesc);
			}
		else if (p->number_of_visits == 1) {
			p->number_of_visits++;
			newicktree+=",";
			Print_MCcoal_Tree(p->rightdesc);
			}
		else if (p->number_of_visits == 2) {  //All descendants have been visited
			p->number_of_visits = 0;
			newicktree+=")";
			if ((p->brlength != -999) || ((p->isroot) && (p->leftdesc->brlength != -999))) {
				newicktree+=(":"+to_string((getnodeage(p)*this->locusrate)));
				}
			if (p->popsize != -999) {
				newicktree+=(" #"+to_string((p->popsize*this->heredityscalar*this->locusrate)));
				}
			if (!p->isroot) {
				Print_MCcoal_Tree(p->ancestor);
				}
			else {
				newicktree+=";";
				//return(void);
				}
			}
		}
	}

	
void Tree::Print_Tree(Node* p)
	{
	if (p->tip) {
		newicktree+=p->species;
		if (p->brlength != -999) {
			newicktree+=(":"+to_string(p->brlength));
			}
		if (p->popsize != -999) {
			newicktree+=(" #"+to_string(p->popsize));
			}
		Print_Tree(p->ancestor);
		}
	else {
		if (p->number_of_visits == 0) {
			p->number_of_visits++;
			newicktree+="(";
			Print_Tree(p->leftdesc);
			}
		else if (p->number_of_visits == 1) {
			p->number_of_visits++;
			newicktree+=",";
			Print_Tree(p->rightdesc);
			}
		else if (p->number_of_visits == 2) {  //All descendants have been visited
			p->number_of_visits = 0;
			newicktree+=")";
			if (p->brlength != -999) {
				newicktree+=(":"+to_string(p->brlength));
				}
			if (p->popsize != -999) {
				newicktree+=(" #"+to_string(p->popsize));
				}
			if (!p->isroot) {
				Print_Tree(p->ancestor);
				}
			else {
				newicktree+=";";
				//return(void);
				}
			}
		}
	}

void Tree::parseTree(string a_tree_string)
{
	tree_string = a_tree_string;
	this->cursor=0;
//cout << "tree_string: " << tree_string << endl;
	while (tree_string[this->cursor] != '(') {
		this->cursor++; 
		}
	this->cursor++;                       //   -> Initialize and recursive call using the basal node
	Node* newnode = new Node;           // Add node
	root = newnode;
	root->isroot=true;
	nextnode(root);
}

void Tree::nextnode(Node*& actual_node)
{
	string temp;
	if (tree_string[this->cursor] == '(') {   //the left node is a node
		this->cursor++;
		Node* newnode = new Node;           // Add node
		actual_node->leftdesc = newnode;
		newnode->ancestor = actual_node;
		nextnode(newnode);      // Go to next node
		}
	if (!strchr("():;#,",tree_string[this->cursor])) {        //the left node is a leaf
		temp="";
		while (!strchr("():;#,[]",tree_string[this->cursor])) {
			temp+=tree_string[this->cursor];
			this->cursor++;
			}
//cout << "Species: " << temp << endl;
		Node* newnode = new Node;             // Add node (tip)
		actual_node->leftdesc = newnode;
		newnode->ancestor = actual_node;
		newnode->tip = true;
		newnode->species = temp;
		if (tree_string[this->cursor] == '[') {             // Population size
			temp="";
			this->cursor++;
			while (!strchr("=",tree_string[this->cursor])) {
				this->cursor++;
				}
			this->cursor++;
			while (!strchr("]",tree_string[this->cursor])) {
				temp+=tree_string[this->cursor];
				this->cursor++;
				}
			this->cursor++;
			newnode->popsize = string_to_double(temp);
//cout << "  popsize: " << temp << endl;
			}
		if (tree_string[this->cursor] == ':') {             // Branch length
			temp="";
			this->cursor++;
			while (!strchr("():;#,",tree_string[this->cursor])) {
				temp+=tree_string[this->cursor];
				this->cursor++;
				}
			newnode->brlength = string_to_double(temp);
//cout << "  brlength: " << temp << endl;
			}
		}
	if (tree_string[this->cursor] == ',') { // Look for things related to the right node
		this->cursor++;
		if (tree_string[this->cursor] == '(') { //the right node is a node ... recursive call
			this->cursor++;
			Node* newnode = new Node;             // Add node
			actual_node->rightdesc = newnode;
			newnode->ancestor = actual_node;
			nextnode(newnode);
			}
		if (!strchr("():;#,",tree_string[this->cursor])) { //the right node is a leaf
			temp="";
			while (!strchr("():;#,[]",tree_string[this->cursor])) {
				temp+=tree_string[this->cursor];
				this->cursor++;
				}
//cout << "Species: " << temp << endl;
			Node* newnode = new Node;              // Add node
			actual_node->rightdesc = newnode;
			newnode->ancestor = actual_node;
			newnode->tip = true;
			newnode->species = temp;
			if (tree_string[this->cursor] == '[') {             // Population size
				temp="";
				this->cursor++;
				while (!strchr("=",tree_string[this->cursor])) {
					this->cursor++;
					}
				this->cursor++;
				while (!strchr("]",tree_string[this->cursor])) {
					temp+=tree_string[this->cursor];
					this->cursor++;
					}
				this->cursor++;
				newnode->popsize = string_to_double(temp);
//cout << "  popsize(right): " << temp << endl;
				}
			if (tree_string[this->cursor] == ':') {             // Branch length
				temp="";
				this->cursor++;
				while (!strchr("():;#,",tree_string[this->cursor])) {
					temp+=tree_string[this->cursor];
					this->cursor++;
					}
				newnode->brlength = string_to_double(temp);
//cout << "  brlength(right): " << temp << endl;
				}
			}
		}
	if (tree_string[this->cursor] == ')') {              // End of node
		this->cursor++;     	
		if (tree_string[this->cursor] == '[') {             // Population size
			temp="";
			this->cursor++;
			while (!strchr("=",tree_string[this->cursor])) {
				this->cursor++;
				}
			this->cursor++;
			while (!strchr("]",tree_string[this->cursor])) {
				temp+=tree_string[this->cursor];
				this->cursor++;
				}
			this->cursor++;
			actual_node->popsize = string_to_double(temp);
//cout << "  popsize(node): " << temp << endl;
			}
		if (tree_string[this->cursor] == ':') {             // Branch length
			temp="";
			this->cursor++;
			while (!strchr("():;#,",tree_string[this->cursor])) {
				temp+=tree_string[this->cursor];
				this->cursor++;
				}
			actual_node->brlength = string_to_double(temp);
			}
		}
	else if (tree_string[this->cursor] == ';') {              // End of tree
		}
} 

	
#define _Tree_CLASS_DEFINED_
#endif /* _Tree_CLASS_DEFINED_ */

