//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains functions that handle the computation of the dimensionalities of a system of beta-functions
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _DIMENSIONS_H_INCLUDED 
#define _DIMENSIONS_H_INCLUDED
#include "common.h"
#include "containers.h"

vector<vector<int> > findDims(const vector<container> &numericalBetaFuncs_1, const vector<container> &numericalBetaFuncs_2, vector<int> &dimFactors);
void findSearchDims(vector<vector<int> > &searchDims, const int &size, const vector<int> &curDims);

#endif

