//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains the overhead functions that set the monomial, polynomial and factorized polynomial searching algorithms in motion. 
//It also contains the overhead for those functions, the findInvariants function, which is called from main.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "invariants.h"
#include "common.h"
#include "symbolics.h"
#include "containers.h"
#include "combinatorics.h"
#include "matrix.h"
#include "linear.h"
#include "quadratic.h"
#include "dimensions.h"
#include "filtering.h"


vector<container> findMono(const vector<container> &numericalBetaFuncs, const vector<vector<int> >&dims)
//The algorithm for finding monomial invariants
{
	map<__uint128_t, unsigned long> valMap;									//Map a combination of powers with a coefficient
	vector<vector<MatEl> > M;
	vector<MatEl> MatRow;
	unsigned long aVal = 0;											//Keep track of the used Coefficients
	for (unsigned long i=0; i<numericalBetaFuncs.size(); i++)						//Convert to a matrix
	{
		for (unsigned long j=0; j<numericalBetaFuncs[i].Coefficients.size(); j++)
		{
			vector<Term> tempAddTo;
			tempAddTo.push_back(Term(-1,i));
			__uint128_t pairNum = termCantor(numericalBetaFuncs[i].Powers[j] + tempAddTo);

			map<__uint128_t, unsigned long>::iterator it;
			it = valMap.find(pairNum);
			if (it==valMap.end())
			{
				valMap[pairNum] = aVal;
				M.push_back(MatRow);
				M[aVal].push_back(MatEl(numericalBetaFuncs[i].Coefficients[j],i));
				aVal++;
			}
			else
			{
				unsigned long oldRow = it->second;
				if (M[oldRow].size()!=0 && M[oldRow].back().col==i)
				{
					M[oldRow].back().val+=numericalBetaFuncs[i].Coefficients[j];
					if (M[oldRow].back().val==0){M[oldRow].pop_back();}
				}
				else {M[oldRow].push_back(MatEl(numericalBetaFuncs[i].Coefficients[j], i));}
			}
		}
	}

	for (unsigned long i=0; i<M.size(); i++){rowGCD(M[i]);}

	vector<vector<long> > nullSpace = findNull (M, numericalBetaFuncs.size());				//Determine nullspace

	vector<container> invariants;										//Convert back to invariants
	for (unsigned long i=0; i<nullSpace.size(); i++)
	{
		container newInv;
		vector<Term> temp;
		vector<int> newDims(dims.size());
		for (unsigned long j=0; j<numericalBetaFuncs.size(); j++)
		{
			if (nullSpace[i][j]!=0)
			{
				temp.push_back(Term(nullSpace[i][j], j));
				for (int k=0; k<dims.size(); k++)
				{
					newDims[k] += nullSpace[i][j] * dims[k][j];
				}
			}
		}
		newInv.addTerm(1,temp);
		newInv.dims = newDims;
		invariants.push_back(newInv);
	}

	for (unsigned long i=0; i<invariants.size(); i++)
	{
		vector<int> cDims(dims.size());
		for (unsigned long j=0; j<invariants[i].Powers[0].size(); j++)
		{
			for (unsigned long k=0; k<dims.size(); k++)
			{
				int index = invariants[i].Powers[0][j].var;
				cDims[k] += invariants[i].Powers[0][j].pow * dims[k][index];
			}
		}
		invariants[i].dims = cDims;
	}

	return invariants;
}

void findLinear(const vector<container> &BetaFuncs_1, const vector<container> &BetaFuncs_2, vector<int> searchDims, vector<vector<int> > dims, vector<int> correctionFactors, vector<container> &newInvariants, map<int8_t,string> symbols)
//The algorithm for finding polynomial invariants.
{
	if(REPORT)
	{
		cout << "Linear search at dims: ";
		for (int i=0; i<searchDims.size(); i++)
		{
			cout << i << ": " << searchDims[i] << " ";
		}
		cout << endl << "Doing Combinatorics " << endl;
	}

	vector<vector<Term> > allTerms_1;								//Contains all possible terms for the 1-loop part of the invariant
	for (unsigned int i=1; i<=MAX_TERM; i++){makeTerms(allTerms_1, dims, searchDims, i);}

	vector<vector<Term> > allTerms = allTerms_1;							//Also make a list of all terms together
	vector<vector<Term> > allTerms_2;								//Contains all possible terms for the 2-loop part of the invariant
	if (INCLUDE_TWO_LOOP)
	{
		vector<int> twoLoopSearchDims = searchDims;
		for (int i=0; i<searchDims.size(); i++)
		{
			twoLoopSearchDims[i] += correctionFactors[i];
		}
		for (unsigned int i=1; i<=MAX_TERM; i++){makeTerms(allTerms_2, dims, twoLoopSearchDims, i);}
		allTerms.insert(allTerms.end(), allTerms_2.begin(), allTerms_2.end());			//Add the two-loop terms to the big list
	}


	if(REPORT) {cout << "Creating Matrix" << endl;}
	createLinearMat(BetaFuncs_1, BetaFuncs_2, allTerms_1, allTerms_2);				//Create the Matrix && write to file



	vector<vector<MatEl> > M;
	if (REPORT) {cout << "Reading" << endl;}
	readMat("Matrix.bin", M);									//Read the matrix from the file
	for (unsigned long i=0; i<M.size(); i++){rowGCD(M[i]);}

	unsigned long rows = M.size();
	unsigned long cols = allTerms.size();

	unsigned long totalElements=0;
	for (unsigned long i=0; i<M.size(); i++){totalElements+=M[i].size();}
	if(REPORT)
	{
		cout << "Done creating matrix: " << rows << " x " << cols << endl;
		cout << "Total nonzero elements: " << totalElements << endl;
		cout << "Sparsity: " << (double)totalElements/((double)(rows*cols)) << endl;
		cout << "Finding nullspace " << endl;
	}

	vector<vector<long> > nullSpace = findNull (M, cols); 						//Find the null space

	vector<container> invarVec;
	for (unsigned long i=0; i<nullSpace.size(); i++)
	{
		container newInv(searchDims);
		for (unsigned long j=0; j<allTerms.size(); j++)						//Fill Polynomial part
		{
			if (nullSpace[i][j]!=0)
			{
				newInv.addTerm(nullSpace[i][j], allTerms[j]);
			}
		}
		invarVec.push_back(newInv);
	}

	for (unsigned long i=0; i<invarVec.size(); i++)							//Reduce the invariants by their gcd
	{
		invarVec[i].reduce();
		if (INCLUDE_TWO_LOOP)									//Filter out one-loop invariants
		{
			if (!invarVec[i].checkDims(dims)){newInvariants.push_back(invarVec[i]);}
		}

		else {newInvariants.push_back(invarVec[i]);}
	}


	remove("Matrix.bin");
	if(REPORT) {cout << "Done" << endl << endl;}
}


void findQuadratic (const vector<container> &BetaFuncs, vector<int> searchDims, vector<vector<int> > dims, vector<container> &newInvariants, map<vector<int>, long> &dimDatabase)
//The algorithm for finding factorized polynomial ivnariants
{
	if(REPORT)
	{
		cout << "Quadratic search at dims ";
		for (int i=0; i<searchDims.size(); i++)
		{
			cout << i << ":" << searchDims[i] << " ";
		}
		cout << endl << "Doing Combinatorics " << endl;
	}
	vector<vector<Term> > allTerms;									//Contains all possible terms
	if(REPORT) {cout << "Doing Combinatorics " << endl;}
	for (unsigned int i=1; i<=MAX_TERM; i++){makeTerms(allTerms, dims, searchDims, i);}

	if(REPORT) {cout << "Creating Matrices" << endl;}
	createQuadraticMat(BetaFuncs, allTerms);							//Create linear & quadratic matrices


	for (int i=0; i<BetaFuncs.size(); i++)
	//for (int i=0; i<1; i++)
	{
		if(REPORT) {cout << "Factorize variable " << i << endl;}
		unsigned long cols = allTerms.size();
		vector<vector<polynomial> > smallMat = createSystem(i, cols);				//Find the accompanying system
		vector<long> eigenvalues = getEigenvaluesNew(smallMat, cols);

		map<vector<int>, long> newdims;								//Database for the new dims
		for (int j=0; j<eigenvalues.size(); j++)
		{
			vector<int> newKey = searchDims;
			for (int k=0; k<dims.size(); k++)
			{
				searchDims[k] += eigenvalues[j]*dims[k][i];
			}
			map<vector<int>, long>::iterator it;
			it  = newdims.find(newKey);
			if (it==newdims.end()) {newdims[newKey] = 1;}
			else {(it->second)++;}
		}
		vector<long> newEigenvalues;								//Find all eigenvalues that correspond with new invariants
		map<vector<int>, long>::iterator it;
		for (it=newdims.begin(); it!=newdims.end(); it++)
		{
			vector<int> key = it->first;
			map<vector<int>, long>::iterator databaseIt;
			databaseIt = dimDatabase.find(key);
			if (databaseIt==dimDatabase.end() || (it->second) > (databaseIt->second))
			{
				for (int j=0; j<dims.size(); j++)					//reconstruct the eigenvalue
				{
					if (dims[j][i]!=0)
					{
						newEigenvalues.push_back((it->first)[j]-searchDims[j]/dims[j][i]);
						break;
					}
				}
			}
		}
		for (int j=0; j<newEigenvalues.size(); j++)						//Find the invariants
		{
			vector<vector<MatEl> > EVM;
			readQuadraticMat(i, newEigenvalues[j], EVM);
			for (unsigned long k=0; k<EVM.size(); k++){rowGCD(EVM[k]);}

			vector<vector<long> > nullSpace = findNull (EVM, allTerms.size()); 		//Find the null space

			vector<container> invarVec;
			for (unsigned long k=0; k<nullSpace.size(); k++)
			{
				vector<int> newDims = searchDims;
				for (int l=0; l<dims.size(); l++)
				{
					newDims[l] += newEigenvalues[j]*dims[l][i];
				}
				container newInv(newDims);
				for (unsigned long l=0; l<allTerms.size(); l++)				//Fill Polynomial part
				{
					if (nullSpace[k][l]!=0)
					{
						vector<Term> newPows = allTerms[l];
						long index = 0;						//Add the eigenvalue to the powers
						while(newPows[index].var < i && index<(newPows.size()-1)){index++;}
						if (newPows[index].var == i)
						{
							newPows[index].pow += newEigenvalues[j];
							if (newPows[index].pow==0) {newPows.erase(newPows.begin() + index);}
						}
						else {newPows.insert(newPows.begin() + index, Term(newEigenvalues[j], i));}
						newInv.addTerm(nullSpace[k][l], newPows);
					}
				}
				invarVec.push_back(newInv);
			}

			for (unsigned long k=0; k<invarVec.size(); k++)					//Reduce the invariants by their gcd
			{
				if (invarVec[k].isMono()==false)
				{
					invarVec[k].reduce();
					newInvariants.push_back(invarVec[k]);
				}
			}
		}

		string name = "Matrix";
		ostringstream tempStream;
		tempStream << i;
		name.append(tempStream.str());
		name.append(".bin");									//Remove the quadratic matrix
		remove(name.c_str());
	}

	remove("Matrix.bin");
	if(REPORT) {cout << "Done" << endl << endl;}
}

void findInvariants(vector<BetaFunc> funcs_1, vector<BetaFunc> funcs_2)
//Master function that converts the Beta-functions from symbolic to numeric entities. Calls the above algorithms.
{
	map<int8_t,string> symbols;									//Association between names and numbers
	map<Variable,int8_t> ordering;
	long N=0;
	for (unsigned long i=0; i<funcs_1.size(); i++)
	{
		for (unsigned long j=0; j<(funcs_1[i].par.size)*(funcs_1[i].par.size); j++)		//Of all elements of the matrices
		{
			if (funcs_1[i].par.vars[j].isZero==false && ordering.count(funcs_1[i].par.vars[j])==false)
			{
				symbols[N] = funcs_1[i].par.vars[j].name;				//Create number association
				if (funcs_1[i].par.vars[j].isComplex==true && funcs_1[i].par.vars[j].CC==true){symbols[N].append("c");}
				ordering[funcs_1[i].par.vars[j]] = N;
				N++;
			}
		}
	}

	if (REPORT) {cout << "There are " << N << " parameters" << endl;}
	long factor_1=1;
	for (unsigned long i=0; i<funcs_1.size(); i++)							//Covert to integer equations
	{
		for (unsigned long j=0; j<(funcs_1[i].par.size)*(funcs_1[i].par.size); j++)
		{
			for (unsigned long k=0; k<funcs_1[i].func.terms[j].size(); k++)
			{
				factor_1 = abs(factor_1*funcs_1[i].func.coeffs[j][k].den)/abs(gcd(factor_1,funcs_1[i].func.coeffs[j][k].den));		//Calculate lcm
			}
			if (INCLUDE_TWO_LOOP)
			{
				for (unsigned long k=0; k<funcs_2[i].func.terms[j].size(); k++)
				{
					factor_1 = abs(factor_1*funcs_2[i].func.coeffs[j][k].den)/abs(gcd(factor_1,funcs_2[i].func.coeffs[j][k].den));	//Calculate lcm
				}
			}
		}
	}

	map<Variable,long> checkDoubles;
	vector<container> numericalBetaFuncs_1;								//Store betafuncs in numerical form
	vector<container> numericalBetaFuncs_2;
	for (unsigned long i=0; i<funcs_1.size(); i++)							//Convert the Matrix-form, symbolic BetaFuncs to numerical funcs
	{
		for (unsigned long j=0; j<(funcs_1[i].par.size)*(funcs_1[i].par.size); j++)
		{
			if (funcs_1[i].par.vars[j].isZero==false && checkDoubles.count(funcs_1[i].par.vars[j])==false)
			{
				container temp_1;							//One loop part
				for (unsigned long k=0; k<funcs_1[i].func.terms[j].size(); k++)
				{
					funcs_1[i].func.coeffs[j][k] = funcs_1[i].func.coeffs[j][k]*ir(factor_1);
					vector<Term> tempvec_1;
					for (unsigned long l=0; l<funcs_1[i].func.terms[j][k].size(); l++)
					{
						vector<Term> addToVec;
						addToVec.push_back(Term(1,ordering[funcs_1[i].func.terms[j][k][l]]));
						tempvec_1 = tempvec_1 + addToVec;
					}
					temp_1.addTermSafe(funcs_1[i].func.coeffs[j][k].num, tempvec_1);
				}
				numericalBetaFuncs_1.push_back(temp_1);
				checkDoubles[funcs_1[i].par.vars[j]] = 1;

				if(INCLUDE_TWO_LOOP)
				{
					container temp_2;						//Two loop part
					for (unsigned long k=0; k<funcs_2[i].func.terms[j].size(); k++)
					{
						funcs_2[i].func.coeffs[j][k] = funcs_2[i].func.coeffs[j][k]*ir(factor_1);
						vector<Term> tempvec_2;
						for (unsigned long l=0; l<funcs_2[i].func.terms[j][k].size(); l++)
						{
							vector<Term> addToVec;
							addToVec.push_back(Term(1,ordering[funcs_2[i].func.terms[j][k][l]]));
							tempvec_2 = tempvec_2 + addToVec;
						}
						temp_2.addTermSafe(funcs_2[i].func.coeffs[j][k].num, tempvec_2);
					}
					numericalBetaFuncs_2.push_back(temp_2);
					checkDoubles[funcs_2[i].par.vars[j]] = 1;
				}
			}
		}
	}

	if(REPORT)
	{
		for (unsigned long i=0; i<numericalBetaFuncs_1.size(); i++)
		{
			cout << i << ": Beta(" << symbols[i] << ") = ";
			numericalBetaFuncs_1[i].print(symbols);
			if (INCLUDE_TWO_LOOP)
			{
				cout << "+ 1/(16pi²) ";
				numericalBetaFuncs_2[i].print(symbols);
			}
		}
	}

	vector<vector<int> > dims;									//Generate dimensions
	vector<int> dimFactors;
	dims = findDims(numericalBetaFuncs_1, numericalBetaFuncs_2, dimFactors);

	if (REPORT)
	{
		cout << "Found the following dims: " << endl;
		for (int i=0; i<dims.size(); i++)
		{
			for (int j=0; j<dims[i].size(); j++)
			{
				cout << dims[i][j] << " ";
			}
			cout << endl;
		}
	}
	vector<vector<int> > searchDims;								//Get all dims and couplDims to search
	vector<int> empty;
	findSearchDims(searchDims, dims.size(), empty);

	map<vector<int>, long> dimDatabase;								//Keeps track of how many invariants of a specific dim were found
																//Here we start searching for invariants
	vector<container> uniqueInvariants;								//Store all unique invariants
	vector<container> newInvariants;								//Store all the invariants we find
	vector<container> monomials;									//Store all momomial invariants
	if (INCLUDE_TWO_LOOP==false)									//Search for Monomials first
	{
		cout << "Monomial Searching" << endl;
		monomials = findMono(numericalBetaFuncs_1, dims);
	}

	for (int i=0; i<monomials.size(); i++)
	{
		monomials[i].print(symbols);
	}
	//------------------------------------------Add inverse of monomials for filtering system------------------------------
	/*
	for (unsigned long i=0; i<monomials.size(); i++)
		{
		for (unsigned long j=0; j<monomials[i].Powers[0].size(); j++) {monomials[i].Powers[0][j].pow = -monomials[i].Powers[0][j].pow;}
		monomials[i].dim = -monomials[i].dim;
		monomials[i].couplDim=-monomials[i].couplDim;
		uniqueInvariants.push_back(monomials[i]);
		}
	*/
	//----------------------------------------------------------------------------------------------------------------------

	long dimIndex=0;
	if (REPORT) {cout << "Polynomial Searching" << endl;}
	while(dimIndex < searchDims.size())								//First find all invariants without factorizing
	{
		findLinear(numericalBetaFuncs_1, numericalBetaFuncs_2, searchDims[dimIndex], dims, dimFactors, newInvariants, symbols);
		dimIndex++;
	}

	for (int i=0; i<newInvariants.size(); i++)							//Add all new invariants to the database
	{
		vector<int> key = newInvariants[i].dims;
		map<vector<int>, long>::iterator it;
		it = dimDatabase.find(key);
		if (it!=dimDatabase.end()) {(it->second)++;}						//Add one to the count if the key already existed
		else {dimDatabase[key]=1;}								//Otherwise add it to the database
	}

	if (FACTORIZE==true && INCLUDE_TWO_LOOP==false)
	{
		if (REPORT) {cout << "Factorized Polynomial Searching" << endl;}
		dimIndex=0;
		while(dimIndex < searchDims.size())
		{
			findQuadratic(numericalBetaFuncs_1, searchDims[dimIndex], dims, newInvariants, dimDatabase);
			dimIndex++;
		}
	}

	sort(newInvariants.begin(), newInvariants.end(), compInv);					//Sorting the new invariants accordingly
	for (int i=0; i<monomials.size(); i++) {newInvariants.insert(newInvariants.begin(), monomials[i]);}

	if (INCLUDE_TWO_LOOP==false)
	{
		if (REPORT) {cout << "Filtering" << endl;}
		filterInvariants(uniqueInvariants, newInvariants, symbols, dims.size());
	}
	else
	{
		uniqueInvariants = newInvariants;
	}
	if (REPORT) {cout << "Done" << endl;}
	for (unsigned long i=0; i<uniqueInvariants.size(); i++){uniqueInvariants[i].dimPrint(symbols, dims);}
}
