//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The Functions in this file create and read the large, sparse matrices that are used in the computation. To save memory, the matrices are first written to 
//the hard drive. Afterwards, when the data structures that are required to build the matrices (the map defined in the first line of the create functions) 
//are deleted. the matrices are read using the read functions.
//Matrices involved in regular polynomial searching are referred to as linear (since they are linear systems), while those involved in factorized polynomial 
//searching are referred to as quadratic (since the system is quadratic in the unknowns).
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "matrix.h"
#include "common.h"
#include "../settings.h"
#include "containers.h"

void createLinearMat(const vector<container> &BetaFuncs_1, const vector<container> &BetaFuncs_2, vector<vector<Term> > &allTerms_1, vector<vector<Term> > &allTerms_2)
//Creates the matrix used in regular polynomial searching
{
	map<__uint128_t, unsigned long> valMap;									//Map a combination of powers with a coefficient
	ofstream fout("Matrix.bin", ios::binary);
	unsigned long aVal = 0;											//Keep track of the used Coefficients

	for (unsigned int i=0; i<allTerms_1.size(); i++)							//for all possible terms
	{
		for (unsigned int j=0; j<allTerms_1[i].size(); j++)						//Loot at each power in those terms
		{
			int currentVar = allTerms_1[i][j].var;
			for (unsigned int k=0; k<BetaFuncs_1[currentVar].Coefficients.size(); k++)		//Look at the terms inside the corresponding beta function
			{
				vector<Term> tempAddTo; 
				tempAddTo.push_back(Term(-1,currentVar));
				__uint128_t pairNum = termCantor(allTerms_1[i] + BetaFuncs_1[currentVar].Powers[k] + tempAddTo);

				long coeff = static_cast<long>(allTerms_1[i][j].pow)*BetaFuncs_1[currentVar].Coefficients[k];
				fout.write(reinterpret_cast<char*>(&coeff), 8);

				map<__uint128_t, unsigned long>::iterator it;
				it = valMap.find(pairNum);
				if (it==valMap.end())								//Check if we've seen the combination of powers in that derivative before
				{
					valMap[pairNum] = aVal;
					fout.write(reinterpret_cast<char*>(&aVal), 8);
					aVal++;
				}
				else 
				{
					unsigned long row = it->second;
					fout.write(reinterpret_cast<char*>(&row), 8);
				}
				fout.write(reinterpret_cast<char*>(&i), 4);
			}
		}
	}

//From here, the two-loop beta-functions are also included if the used has indicated that they should be. Loop over both sets of terms again 
//since the two-loop effects originate from both the one-loop and two-loop terms in the invariant.
	if(INCLUDE_TWO_LOOP)
	{
		for (unsigned int i=0; i<allTerms_1.size(); i++)				
		{
			for (unsigned int j=0; j<allTerms_1[i].size(); j++)			
			{
				int currentVar = allTerms_1[i][j].var;
				for (unsigned int k=0; k<BetaFuncs_2[currentVar].Coefficients.size(); k++)	
				{
					vector<Term> tempAddTo; 
					tempAddTo.push_back(Term(-1,currentVar));
					__uint128_t pairNum = termCantor(allTerms_1[i] + BetaFuncs_2[currentVar].Powers[k] + tempAddTo);

					long coeff = static_cast<long>(allTerms_1[i][j].pow)*BetaFuncs_2[currentVar].Coefficients[k];
					fout.write(reinterpret_cast<char*>(&coeff), 8);

					map<__uint128_t, unsigned long>::iterator it;
					it = valMap.find(pairNum);
					if (it==valMap.end())					
					{
						valMap[pairNum] = aVal;
						fout.write(reinterpret_cast<char*>(&aVal), 8);
						aVal++;
					}
					else 
					{
						unsigned long row = it->second;
						fout.write(reinterpret_cast<char*>(&row), 8);
					}
					fout.write(reinterpret_cast<char*>(&i), 4);
				}
			}
		} 

		for (unsigned int i=0; i<allTerms_2.size(); i++)						//for all possible terms
		{
			for (unsigned int j=0; j<allTerms_2[i].size(); j++)					//Loot at each power in those terms
			{
				int currentVar = allTerms_2[i][j].var;
				for (unsigned int k=0; k<BetaFuncs_1[currentVar].Coefficients.size(); k++)	//Look at the terms inside the corresponding beta function
				{
					vector<Term> tempAddTo; 
					tempAddTo.push_back(Term(-1,currentVar));
					__uint128_t pairNum = termCantor(allTerms_2[i] + BetaFuncs_1[currentVar].Powers[k] + tempAddTo);

					long coeff = static_cast<long>(allTerms_2[i][j].pow)*BetaFuncs_1[currentVar].Coefficients[k];
					fout.write(reinterpret_cast<char*>(&coeff), 8);

					map<__uint128_t, unsigned long>::iterator it;
					it = valMap.find(pairNum);
					if (it==valMap.end())							//Check if we've seen the combination of powers in that derivative before
					{
						valMap[pairNum] = aVal;
						fout.write(reinterpret_cast<char*>(&aVal), 8);
						aVal++;
					}
					else 
					{
						unsigned long row = it->second;
						fout.write(reinterpret_cast<char*>(&row), 8);
					}
					long col = i+allTerms_1.size();
					fout.write(reinterpret_cast<char*>(&col), 4);
				}
			}
		}
	}
	fout.close();
}

void createQuadraticMat(const vector<container> &BetaFuncs, vector<vector<Term> > &allTerms)
//This function creates the matrices required for factorized polynomial searching. It includes the same procedure as the 
//linear function, but then continues to build a matrix for the factorization of every parameter
{
	map<__uint128_t, unsigned long> valMap;									//Because the map is universal, we have to create all matrices at once
	ofstream fout("Matrix.bin", ios::binary);
	unsigned long aVal = 0;											//Keep track of the used Coefficients

	for (unsigned int i=0; i<allTerms.size(); i++)								//for all possible terms
	{
		for (unsigned int j=0; j<allTerms[i].size(); j++)						//Loot at each power in those terms
		{
			int currentVar = allTerms[i][j].var;
			for (unsigned int k=0; k<BetaFuncs[currentVar].Coefficients.size(); k++)		//Look at the terms inside the corresponding beta function
			{
				vector<Term> tempAddTo; 
				tempAddTo.push_back(Term(-1,currentVar));
				__uint128_t pairNum = termCantor(allTerms[i] + BetaFuncs[currentVar].Powers[k] + tempAddTo);

				long coeff = static_cast<long>(allTerms[i][j].pow)*BetaFuncs[currentVar].Coefficients[k];
				fout.write(reinterpret_cast<char*>(&coeff), 8);

				map<__uint128_t, unsigned long>::iterator it;
				it = valMap.find(pairNum);
				if (it==valMap.end())								//Check if we've seen the combination of powers in that derivative before
				{
					valMap[pairNum] = aVal;
					fout.write(reinterpret_cast<char*>(&aVal), 8);
					aVal++;
				}
				else 
				{
					unsigned long row = it->second;
					fout.write(reinterpret_cast<char*>(&row), 8);
				}
				fout.write(reinterpret_cast<char*>(&i), 4);
			}
		}
	}

//From here, the matrices for factorization are built
	fout.close();
	if (FACTORIZE)
	{
		for (unsigned int i=0; i<BetaFuncs.size(); i++)							//Loop over all parameters that can be factorized
		{
			string name = "Matrix";  								//Create file name
			ostringstream tempStream;
			tempStream << i;
			string tempString = tempStream.str();
			name.append(tempString);
			name.append(".bin");
			
			ofstream fout(name, ios::binary);
			for (unsigned int j=0; j<allTerms.size(); j++)						//for all possible terms
			{
				for (unsigned long k=0; k<BetaFuncs[i].Coefficients.size(); k++)		//Look at the terms inside the corresponding beta function
				{
					vector<Term> tempAddTo; 
					tempAddTo.push_back(Term(-1, i));
					__uint128_t pairNum = termCantor(allTerms[j] + BetaFuncs[i].Powers[k] + tempAddTo);

					long coeff = BetaFuncs[i].Coefficients[k];
					fout.write(reinterpret_cast<char*>(&coeff), 8);

					map<__uint128_t, unsigned long>::iterator it;
					it = valMap.find(pairNum);

					if (it==valMap.end())							//Check if we've seen the combination of powers in that derivative before
					{
						valMap[pairNum] = aVal;
						fout.write(reinterpret_cast<char*>(&aVal), 8);
						aVal++;
					}
					else 
					{
						unsigned long row = it->second;
						fout.write(reinterpret_cast<char*>(&row), 8);
					}	
					fout.write(reinterpret_cast<char*>(&j), 4); 
				}
			}
			fout.close();
		}
	}
}

void readMat(const string &name, vector<vector<MatEl> > &M)
//Read a linear matrix from the hard drive.
{
	ifstream fin(name, ios::binary);
	vector<MatEl> MatRow;
	long coeff;
	unsigned long row;
	unsigned int col;
	while(true)
	{
		fin.read(reinterpret_cast<char*>(&coeff), 8);
		if(fin.eof()){break;}
		fin.read(reinterpret_cast<char*>(&row), 8);
		fin.read(reinterpret_cast<char*>(&col), 4);
		while(M.size()<row){M.push_back(MatRow);}
		if (M.size()==row)
		{
			M.push_back(MatRow);
			M.back().push_back(MatEl(coeff, col));
		}
		else
		{
			if (M[row].size()!=0 && M[row].back().col==col)
			{
				M[row].back().val+=coeff;
				if (M[row].back().val==0){M[row].pop_back();}
			}
			else {M[row].push_back(MatEl(coeff, col));}
		}
	}
	fin.close();
}


void readQuadraticMat(const long factorizedPow, const long eigenvalue, vector<vector<MatEl> > &EVM)
{
//Read a quadratic matrix from the hard drive
	vector<MatEl> MatRow;
	vector<vector<MatEl> > M;
	long coeff;
	unsigned long row;
	unsigned int col;

	ifstream fin("Matrix.bin", ios::binary);
	while(true)
	{
		fin.read(reinterpret_cast<char*>(&coeff), 8);
		if(fin.eof()){break;}
		fin.read(reinterpret_cast<char*>(&row), 8);
		fin.read(reinterpret_cast<char*>(&col), 4);
		while(M.size()<row){M.push_back(MatRow);}
		if (M.size()==row)
		{
			M.push_back(MatRow);
			M.back().push_back(MatEl(coeff, col));
		}
		else
		{
			if (M[row].size()!=0 && M[row].back().col==col)
			{
				M[row].back().val+=coeff;
				if (M[row].back().val==0){M[row].pop_back();}
			}
			else {M[row].push_back(MatEl(coeff, col));}
		}
	}
	fin.close();
	
	string name = "Matrix";  										//Create name
	ostringstream tempStream;
	tempStream << factorizedPow;
	string tempString = tempStream.str();
	name.append(tempString);
	name.append(".bin");
			
	ifstream fintwo(name, ios::binary);
	while(true)
	{
		fintwo.read(reinterpret_cast<char*>(&coeff), 8);
		if(fintwo.eof()){break;}
		fintwo.read(reinterpret_cast<char*>(&row), 8);
		fintwo.read(reinterpret_cast<char*>(&col), 4);
		while(M.size()<row){M.push_back(MatRow);}
		if (M.size()==row)
		{
			M.push_back(MatRow);
			M.back().push_back(MatEl(eigenvalue*coeff, col));						//Add the eigenvalue matrix
		}
		else 		
		{
			if (M[row].size()!=0)
			{
				long index = 0;
				while(M[row][index].col<col && (index < M[row].size())){index++;}
				if (index == M[row].size())
				{
					M[row].push_back(MatEl(eigenvalue*coeff, col));
				}
				else if (M[row][index].col==col)
				{
					M[row][index].val += eigenvalue*coeff;
					if (M[row][index].val==0){M[row].erase(M[row].begin() + index);}
				}
				else 
				{
					M[row].insert(M[row].begin() + index, MatEl(eigenvalue*coeff, col));
				}
			}
			else {M[row].push_back(MatEl(coeff, col));}			
		}
	}
	fintwo.close();

	while(M.size()!=0)
	{
		if (M.back().size()!=0)
		{
			EVM.push_back(M.back());
		}
		M.pop_back();
	}
	//printMat(M);
}
		
