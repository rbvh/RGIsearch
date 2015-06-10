//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains functions that handle the computation of the dimensionalities of a system of beta-functions
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "dimensions.h"
#include "common.h"
#include "../settings.h"
#include "containers.h"
#include "linear.h"

vector<vector<int> > findDims(const vector<container> &numericalBetaFuncs_1, const vector<container> &numericalBetaFuncs_2, vector<int> &dimFactors)
//Use linear system solving to find the dimensionalities of a system of beta-functions
{
	vector<vector<MatEl> > M;
	vector<MatEl> emptyRow;
	for (int i=0; i<numericalBetaFuncs_1.size(); i++)
	{
		for (int j=0; j<numericalBetaFuncs_1[i].Powers.size(); j++)
		{
			M.push_back(emptyRow);
			for (int k=0; k<numericalBetaFuncs_1[i].Powers[j].size(); k++)							//Add all regular terms
			{
				M.back().push_back(MatEl(numericalBetaFuncs_1[i].Powers[j][k].pow, numericalBetaFuncs_1[i].Powers[j][k].var));
			}
			int index=0;													
			while(true)
			{		
				if (index == M.back().size())
				{
					M.back().push_back(MatEl(-1, i));
					break;
				}
				else if (M.back()[index].col == i)
				{
					M.back()[index].val--;
					if (M.back()[index].val==0){M.back().erase(M.back().begin() + index);}
					break;
				}
				else if (M.back()[index].col > i)
				{
					M.back().insert(M.back().begin() + index, MatEl(-1, i));
					break;
				}
				index++;
			}
			M.back().push_back(MatEl(-1, numericalBetaFuncs_1.size()));							//Add the term for the dimensional constant
		}
	}

	if (INCLUDE_TWO_LOOP)
	{	
		for (int i=0; i<numericalBetaFuncs_2.size(); i++)
		{
			for (int j=0; j<numericalBetaFuncs_2[i].Powers.size(); j++)
			{
				M.push_back(emptyRow);
				for (int k=0; k<numericalBetaFuncs_2[i].Powers[j].size(); k++)						//Add all regular terms
				{
					M.back().push_back(MatEl(numericalBetaFuncs_2[i].Powers[j][k].pow, numericalBetaFuncs_2[i].Powers[j][k].var));
				}
				int index=0;												//Add the -1 for the dim of the betaFunc variable
				while(true)
				{		
					if (index == M.back().size())
					{
						M.back().push_back(MatEl(-1, i));
						break;
					}
					else if (M.back()[index].col == i)
					{
						M.back()[index].val--;
						if (M.back()[index].val==0){M.back().erase(M.back().begin() + index);}
						break;
					}
					else if (M.back()[index].col > i)
					{
						M.back().insert(M.back().begin() + index, MatEl(-1, i));
						break;
					}
					index++;
				}
				M.back().push_back(MatEl(-1, numericalBetaFuncs_2.size()+1));						//Add the term for the dimensional constant
			}
		}
	}
	vector<vector<long> > nullSpace;

	if (INCLUDE_TWO_LOOP) 
	{
		nullSpace = findNull(M, numericalBetaFuncs_1.size() + 2);
		for (int i=0; i<nullSpace.size(); i++)
		{
			dimFactors.push_back(nullSpace[i][nullSpace[i].size()-1] - nullSpace[i][nullSpace[i].size()-2]);
		}
	}
	else {nullSpace = findNull(M, numericalBetaFuncs_1.size()+1);}

	vector<vector<int> > dims(nullSpace.size());
	for (int i=0; i<nullSpace.size(); i++)
	{
		for (int j=0; j<numericalBetaFuncs_1.size(); j++)
		{
			dims[i].push_back(nullSpace[i][j]);	
		}
	}
	return dims;
}

void findSearchDims(vector<vector<int> > &searchDims, const int &size, const vector<int> &curDims)
//Determines the dimensions to be searched by the rest of the algorithm
{
	if (curDims.size() == size)
	{
		searchDims.push_back(curDims);
		return;
	}

	for (int i= -DIM_SEARCH_PARAMETER; i<= DIM_SEARCH_PARAMETER; i++)
	{
		vector<int> newDims = curDims;
		newDims.push_back(i);
		findSearchDims(searchDims, size, newDims);
	}
}
