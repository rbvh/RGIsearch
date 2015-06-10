//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains the functions that handle the filtering of newly found invariants. 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _FILTERING_H_INCLUDED 
#define _FILTERING_H_INCLUDED 
#include "common.h"
#include "containers.h"
void newMakeInvariants(const vector<container> &uniqueInvariants, const long &prevInvNum, const container &prevProd, vector<container> &prodInvariants, const vector<int> &dims, const long &allowedHNP, map<int8_t,string> symbols);
bool compInv(const container &one, const container &two);
bool isIndep(vector<container> prevInv, container newInv, map<int8_t,string> symbols);
void filterInvariants(vector<container> &uniqueInvariants, vector<container> &newInvariants, map<int8_t,string> symbols, int dimSize);

#endif


