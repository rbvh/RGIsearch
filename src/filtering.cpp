//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains the functions that handle the filtering of newly found invariants. 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "filtering.h"
#include "common.h"
#include "../settings.h"
#include "containers.h"
#include "linear.h"

void newMakeInvariants(const vector<container> &uniqueInvariants, const long &prevInvNum, const container &prevProd, vector<container> &prodInvariants, const vector<int> &dims, const long &allowedHNP, map<int8_t,string> symbols)
//Recursive function that multiplies invariants to find a set that has the same dimensionalities and characteristics as the candidate invariant
{
	for (int i=prevInvNum; i<uniqueInvariants.size(); i++)
	{
		container curProd = prevProd*uniqueInvariants[i];
		unsigned long newHNP=curProd.HNP();

		if (newHNP<=allowedHNP)
		{
			if (curProd.dims==dims)
			{
				prodInvariants.push_back(curProd);
			}
			bool dimCheck=true;
			for (int j=0; j<dims.size(); j++)
			{
				if (abs(curProd.dims[j]) >= abs(dims[j]) + FILTER_THRESHOLD)
				{
					dimCheck = false;
					break;
				}
			}
			if (dimCheck)
			{
				newMakeInvariants(uniqueInvariants, i, curProd, prodInvariants, dims, allowedHNP, symbols);
			}
		}
	}
}


bool compInv(const container &one, const container &two)
//Comparison function for the amount of terms in an invariant
{
	unsigned long oneHighest = one.HNP();
	unsigned long twoHighest = two.HNP();
	if (oneHighest!=twoHighest) {return (oneHighest<twoHighest);}
	else {return one.Powers.size() < two.Powers.size();}
}


bool isIndep(vector<container> prevInv, container newInv, map<int8_t,string> symbols)
//Determine if an invariant is independent of another set of invariants
{
	map<__uint128_t, unsigned long> valMap;
	vector<vector<MatEl> > invMat;
	vector<MatEl> temp;
	unsigned long aVal = 0;
	for (unsigned long i=0; i<prevInv.size(); i++)							//Add the unique invariants
	{
		for(unsigned long j=0; j<prevInv[i].Coefficients.size(); j++)
		{
			__uint128_t pairNum = termCantor(prevInv[i].Powers[j]);
			map<__uint128_t, unsigned long>::iterator it;
			it = valMap.find(pairNum);
			unsigned long row;
			if (it==valMap.end())
			{
				valMap[pairNum] = aVal;
				invMat.push_back(temp);
				row = aVal;
				aVal++;
			}
			else
			{
				row = it->second;
			}
			invMat[row].push_back(MatEl(prevInv[i].Coefficients[j], i));
		}
	}

	for (unsigned long i=0; i<newInv.Coefficients.size(); i++)					//Add the new invariant
	{
		__uint128_t pairNum = termCantor(newInv.Powers[i]);
		map<__uint128_t, unsigned long>::iterator it;
		it = valMap.find(pairNum);
		if (it==valMap.end())									//Some term exists in the new invariants that wasnt anywhere in the old ones
		{
			return true;
		}
		else
		{
			invMat[it->second].push_back(MatEl(newInv.Coefficients[i], prevInv.size()));
		}
	}

	for (unsigned long i=0; i<newInv.Coefficients.size(); i++)
	{
		__uint128_t pairNum = termCantor(newInv.Powers[i]);
		if (valMap.find(pairNum)==valMap.end()){return true;}
	}


	unsigned long s = findNulldim(invMat, prevInv.size()+1);
	if (s==0){return true;}
	else {return false;}
}

//---------------------------------------------Find set of invariants to remove----------------------------------
void filterInvariants(vector<container> &uniqueInvariants, vector<container> &newInvariants, map<int8_t,string> symbols, int dimSize)
//Use the above functions to filter a set of invariants
{
	for (unsigned long i=0; i<newInvariants.size(); i++)
	{
		vector<container> compareInvariants;
		container constant;
		constant.makeConstant(dimSize);

		newMakeInvariants(uniqueInvariants, 0, constant, compareInvariants, newInvariants[i].dims, newInvariants[i].HNP(), symbols);

		map<container, unsigned long> checkDoubles;
		if (isIndep(compareInvariants, newInvariants[i], symbols))
		{
			uniqueInvariants.push_back(newInvariants[i]);
		}
	}
}
