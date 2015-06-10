//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The Functions in this file create and read the large, sparse matrices that are used in the computation. To save memory, the matrices are first written to 
//the hard drive. Afterwards, when the data structures that are required to build the matrices (the map defined in the first line of the create functions) 
//are deleted. the matrices are read using the read functions.
//Matrices involved in regular polynomial searching are referred to as linear (since they are linear systems), while those involved in factorized polynomial 
//searching are referred to as quadratic (since the system is quadratic in the unknowns).
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _MATRIX_H_INCLUDED 
#define _MATRIX_H_INCLUDED
#include "common.h"
#include "containers.h"
void createLinearMat(const vector<container> &BetaFuncs_1, const vector<container> &BetaFuncs_2, vector<vector<Term> > &allTerms_1, vector<vector<Term> > &allTerms_2);
void createQuadraticMat(const vector<container> &BetaFuncs, vector<vector<Term> > &allTerms);
void readMat(const string &name, vector<vector<MatEl> > &M);
void readQuadraticMat(const long factorizedPow, const long eigenvalue, vector<vector<MatEl> > &EVM);
#endif
