//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains the functions that handle the solving of the linear systems in regular polynomial searching
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _LINEAR_H_INCLUDED 
#define _LINEAR_H_INCLUDED
#include "common.h"
#include "containers.h"

void elimOneRow(vector<vector<MatEl> > &M, vector<vector<unsigned long> > &colM, queue<unsigned long> &elimCandidates, vector<Equation> &N);
void elimMoreRow(vector<vector<MatEl> > &M, vector<MatEl> pivRow, vector<vector<unsigned long> > &rowC, vector<vector<unsigned long> > &colC, vector<vector<unsigned long> > &colM, unsigned long elR, unsigned long pivC, unsigned long pivCoord);
void OneCounts(vector<vector<MatEl> > &M, unsigned long cols, vector<Equation> &N);
void Markowitz(vector<vector<MatEl> > &M, unsigned long cols, vector<Equation> &N);
vector<vector<long> > findNull(vector<vector<MatEl> > &M, unsigned long cols);
unsigned long findNulldim(vector<vector<MatEl> > &M, unsigned long cols);
#endif

