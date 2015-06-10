//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The functions in this file handle the procedure of generating all possible monomial terms that are used in the polynomial searching procedure
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "combinatorics.h"
#include "common.h"
#include "containers.h"
void makeParams(vector<vector<long> > &params, const vector<long> &prev, const unsigned long &length, const long &N)
//This function generates a list of vectors that contain all possible combinations of parameters of length &length. 
{
	if (prev.size()==length){params.push_back(prev);}
	else
	{
		if (prev.size()!=0)
		{
			for (long i=prev[prev.size()-1]+1; i<N; i++)
			{
				vector<long> temp = prev;
				temp.push_back(i);
				makeParams(params,temp,length, N);
			}
		}
		else
		{
			for (long i=0; i<N; i++)
			{
				vector<long> temp = prev;
				temp.push_back(i);
				makeParams(params,temp,length, N);
			}
		}
	}
}

void makePows(vector<vector<long> > &pows, const vector<long> &prev, const unsigned long &length, const long &maxPow)
//This function generates a list of vectors that contain all possible powers p_i of parameters where -maxPow<=p_i<=maxPow
{
	if (prev.size()==length){pows.push_back(prev);}
	else
	{
		for (long i=-maxPow; i<=maxPow; i++)
		{
			if (i!=0)
			{
				vector<long> temp = prev;
				temp.push_back(i);
				makePows(pows, temp,length,maxPow);
			}
		}
	}
}

void makeTerms(vector<vector<Term> > &allTerms, const vector<vector<int> > &dims, const vector<int> &searchDims, const int &maxTerms)
//This function combines the results of the above functions to generate all possible monomial terms. It then checks if they have the correct dimensions.
{
	vector<vector<long> > powers;
	vector<vector<long> > params;

	long maxPower = 3;
	for (int i=0; i<searchDims.size(); i++)
	{
		if (abs(searchDims[i])>maxPower) {maxPower = abs(searchDims[i]);}
	}

	vector<long> empty;
	makePows(powers,empty,maxTerms,maxPower);
	makeParams(params,empty,maxTerms,dims[0].size());
	vector<vector<Term> > result;
	for (unsigned long i=0; i<params.size(); i++)
	{
		for (unsigned long j=0; j<powers.size(); j++)
		{
			vector<int> termDims(dims.size());
			for (long k=0; k<maxTerms; k++)
			{
				for (unsigned long l=0; l<dims.size(); l++)
				{
					termDims[l] += dims[l][params[i][k]]*powers[j][k];
				}
			}
			if (termDims == searchDims)
			{
				vector<Term> temp;
				for (long k=0; k<maxTerms; k++){temp.push_back(Term(powers[j][k], params[i][k]));}
				allTerms.push_back(temp);
			}
		}
	}
}
