//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The functions in this file handle the quadratic solving, e.g. generalized eigenvalue computation for factorized polynomial searching.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _QUADRATIC_H_INCLUDED 
#define _QUADRATIC_H_INCLUDED
#include "common.h"
#include "containers.h"

void reduce(vector<polyMatEl> &row);
unsigned long evaluateRow(const vector<vector<MatEl> > &M, const vector<vector<MatEl> > &N, unsigned long row);
void combinePrimes(vector<long> &primePowers, vector<long> &solutions, long coeff, long curSolution);
vector<long> solvePoly(polynomial &p);
vector<vector<polynomial> > createSystem(const int &factorizedPow, unsigned long &cols);
vector<long> getEigenvaluesNew(vector<vector<polynomial> > &smallMat, unsigned long cols);

#endif
