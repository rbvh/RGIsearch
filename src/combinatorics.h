//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The functions in this file handle the procedure of generating all possible monomial terms that are used in the polynomial searching procedure
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _COMBINATORICS_H_INCLUDED 
#define _COMBINATORICS_H_INCLUDED
#include "common.h"
#include "containers.h"
void makeParams(vector<vector<long> > &params, const vector<long> &prev, const unsigned long &length, const long &N);
void makePows(vector<vector<long> > &pows, const vector<long> &prev, const unsigned long &length, const long &maxPow);
void makeTerms(vector<vector<Term> > &allTerms, const vector<vector<int> > &dims, const vector<int> &searchDims, const int &maxTerms);
#endif
