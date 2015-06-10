//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains the overhead functions that set the monomial, polynomial and factorized polynomial searching algorithms in motion. 
//It also contains the overhead for those functions, the findInvariants function, which is called from main.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _INVARIANTS_H_INCLUDED 
#define _INVARIANTS_H_INCLUDED 

#include "common.h"
#include "../settings.h"
#include "symbolics.h"
#include "containers.h"

vector<container> findMono(const vector<container> &numericalBetaFuncs, const vector<vector<int> >&dims);
void findLinear(const vector<container> &BetaFuncs_1, const vector<container> &BetaFuncs_2, vector<int> searchDims, vector<vector<int> > dims, vector<int> correctionFactors, vector<container> &newInvariants, map<int8_t,string> symbols);
void findQuadratic (const vector<container> &BetaFuncs, vector<int> searchDims, vector<vector<int> > dims, vector<container> &newInvariants, map<vector<int>, long> &dimDatabase);
void findInvariants(vector<BetaFunc> funcs_1, vector<BetaFunc> funcs_2);

#endif

