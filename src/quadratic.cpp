//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The functions in this file handle the quadratic solving, e.g. generalized eigenvalue computation for factorized polynomial searching.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "quadratic.h"
#include "common.h"
#include "containers.h"
#include "matrix.h"
void reduce(vector<polyMatEl> &row)
//Perform several simplifications on a matrix row of polynomials.
{
	queue<long> nums;										//First reduce all numbers by their gcd
	for (int i=0; i<row.size(); i++)
	{
		for (int j=0; j<row[i].poly.terms.size(); j++)
		{
			nums.push(row[i].poly.terms[j].coeff);
		}
	}

	long div = 1;
	if (nums.size()==0){return;}
	else if (nums.size()==1) {div = nums.front();}
	else 
	{
		long first = nums.front();
		nums.pop();
		long second = nums.front();
		nums.pop();
		div = gcd(first, second);
		while(nums.size()!=0)
		{
			div = gcd(div, nums.front());
			nums.pop();
		}
	}

	for (int i=0; i<row.size(); i++)
	{
		for (int j=0; j<row[i].poly.terms.size(); j++)
		{
			row[i].poly.terms[j].coeff = row[i].poly.terms[j].coeff/div;
		}
	}
	unsigned long degree = 2000000000;	
	for (int i=0; i<row.size(); i++)
	{
		for (int j=0; j<row[i].poly.terms.size(); j++)
		{
			if (row[i].poly.terms[j].pow < degree){degree = row[i].poly.terms[j].pow;}
		}
	}
	for (int i=0; i<row.size(); i++)
	{
		for (int j=0; j<row[i].poly.terms.size(); j++)
		{
			row[i].poly.terms[j].pow-=degree;
		}
	}
}

unsigned long evaluateRow(const vector<vector<MatEl> > &M, const vector<vector<MatEl> > &N, unsigned long row)
//Helper function for Gaussian elimination
{
	unsigned long result=0;
	if (row < M.size()) {result += M[row].size();}
	result += N[row].size();
	return result;
}

void combinePrimes(vector<long> &primePowers, vector<long> &solutions, long coeff, long curSolution)
//Find the prime factorization of an integer. Used to solve a polynomial equation
{
	if (coeff==primePowers.size())
	{
		solutions.push_back(curSolution);
		solutions.push_back(-curSolution);
		return;
	}
	for (int i=0; i<=primePowers[coeff]; i++)
	{
		combinePrimes(primePowers, solutions, coeff+1, curSolution*pow(primes[coeff], i));
	}
}

vector<long> solvePoly(polynomial &p)
//Solve a polynomial equation.
{
	long factor = abs(p.terms.back().coeff);
	vector<long> primePowers(primes.size());
	int primeCoeff=0;
	while(factor>1)											//Find prime decomposition of the lowest-degree factor
	{
		while(factor%primes[primeCoeff]!=0)
		{
			primeCoeff++;
			if (primeCoeff==primes.size()){cout << "Error: list of primes too short. Tried to decompose " << factor << endl; exit(1);}
		}  
		primePowers[primeCoeff]++;
		factor /= primes[primeCoeff];
	}
	vector<long> candidateSolutions;
	vector<long> solutions;
	combinePrimes(primePowers, candidateSolutions, 0, 1);
	for (int i=0; i<candidateSolutions.size(); i++)
	{
		if (p.evaluate(candidateSolutions[i])==0) {solutions.push_back(candidateSolutions[i]);}
	}
	return solutions;
}

//----------------------------------------------------------------------------Do Quadratic Solving-----------------------------------------
vector<vector<polynomial> > createSystem(const int &factorizedPow, unsigned long &cols)
//Perform quadratic solving.
{
	vector<vector<MatEl> > M;									//Linear matrix. 
	readMat("Matrix.bin", M);									//Read linear matrix

	vector<vector<MatEl> > N;									//Quadratic matrix
	string name = "Matrix";  									//Read quadratic matrix, create name
	ostringstream tempStream;
	tempStream << factorizedPow;
	string tempString = tempStream.str();
	name.append(tempString);
	name.append(".bin");
	readMat(name, N);
	
	vector<vector<unsigned long> > colM(cols);							//col dom form of the matrices
	vector<vector<unsigned long> > colN(cols);
	queue<unsigned long> elimCandidates;								//Holds rows that are about to be eliminated
	unsigned long eliminationCounter=0;								//Keeps track of the amount of eliminations that are done

	for (unsigned long i=0; i<N.size(); i++)							//Loop over the size of N since N.size()>M.size
	{
		unsigned int totalEl=0; 								//Create colM/ColN and find initial elimination candidates
		if (i<M.size())
		{
			for (unsigned int j=0; j<M[i].size(); j++){colM[M[i][j].col].push_back(i);}
			totalEl +=M[i].size();
		}

		for (unsigned int k=0; k<N[i].size(); k++){colN[N[i][k].col].push_back(i);}
		totalEl += N[i].size();

		if (totalEl==1){elimCandidates.push(i);}						//Add candidates for elimination
	}

	while(elimCandidates.size()!=0)									//Do actual eliminations
	{
		unsigned long curRow = elimCandidates.front();
		unsigned long curCol;
		elimCandidates.pop();
		bool didElimination = false;
		if (curRow<M.size() && M[curRow].size()==1)						//Find the column that is to be eliminated
		{
			curCol = M[curRow].front().col;
			vector<MatEl>().swap(M[curRow]);
			didElimination = true;							
		}
		else if (N[curRow].size()==1)
		{
			curCol = N[curRow].front().col;
			vector<MatEl>().swap(N[curRow]);
			didElimination = true;					
		}

		if (didElimination)									//If we did an elimination (e.g. row was not empty), find and erase all entries of that col
		{
			eliminationCounter++;
			for (int i=0; i<colM[curCol].size(); i++)					//Do M first
			{
				unsigned long elimRow = colM[curCol][i];
				if (elimRow!=curRow)
				{
					unsigned long j=0;
					while(M[elimRow][j].col!=curCol){j++;}
					M[elimRow].erase(M[elimRow].begin() + j);
	
					if (evaluateRow(M, N, elimRow)==1){elimCandidates.push(elimRow);}
				}
			}
		

			for (int i=0; i<colN[curCol].size(); i++)					//Then do N
			{
				unsigned long elimRow = colN[curCol][i];
				if (elimRow!=curRow)
				{
					unsigned long j=0;
					while(N[elimRow][j].col!=curCol){j++;}
					N[elimRow].erase(N[elimRow].begin() + j);

					if (evaluateRow(M, N, elimRow)==1){elimCandidates.push(elimRow);}
				}
			}		
		}
	}

	vector<vector<polynomial> > smallMat;								//This matrix is not in sparse format
	map<unsigned long, unsigned long> colToVar;							//Transformation to a system without eliminated variables
	unsigned long aVal = 0;										//Powers are the first variables, so start after those

	unsigned long curVar;
	for (unsigned long i=0; i<N.size(); i++)	
	{
		vector<polynomial> newRow(cols-eliminationCounter);
		bool isZero=true;
		if (i<M.size())
		{
			for (int k=0; k<M[i].size(); k++)
			{
				map<unsigned long, unsigned long>::iterator it;
				it = colToVar.find(M[i][k].col);					//Find if the column was already assigned to a variable
				if (it==colToVar.end())							//If it hasn't, add it
				{
					curVar = aVal;
					colToVar[M[i][k].col] = aVal;
					aVal++;
				}
				else {curVar = it->second;}						//Else, use that var
				monomial tempMono(M[i][k].val, 0); 					
				newRow[curVar].addTerm(tempMono);
				isZero=false;
			}
		}

		for (int k=0; k<N[i].size(); k++)
		{
			map<unsigned long, unsigned long>::iterator it;
			it = colToVar.find(N[i][k].col);
			if(it==colToVar.end())
			{
				curVar = aVal;
				colToVar[N[i][k].col] = aVal;
				aVal++;
			}
			else {curVar = it->second;}
			monomial tempMono(N[i][k].val, 1);
			newRow[curVar].addTerm(tempMono);
			isZero=false;
		}


		if(isZero==false)										//Check if the row is not completely zero
		{
			smallMat.push_back(newRow);
		}
	}

	map<unsigned long, unsigned long>::iterator it;
	//for (it=colToVar.begin(); it!=colToVar.end(); it++) { cout << it->first <<  "->" << it->second << endl;} 

	for (int i=0; i<smallMat.size(); i++)									//Filter duplicates 
	{
		for (int j=0; j<i; j++)
		{
			if (isEqual(smallMat[i], smallMat[j]))
			{
				smallMat.erase(smallMat.begin()+i);
				i--;
				break;
			}
		}
	} 

	for (int i=0; i<smallMat.size(); i++)
	{
		reduceNum(smallMat[i]);
		//print(smallMat[i]);
	}

	cols = aVal;												//We change cols to the amount of cols that are not eliminated
	return smallMat;	
}



vector<long> getEigenvaluesNew(vector<vector<polynomial> > &smallMat, unsigned long cols)
//Compute eigenvalues
	{
	vector<polyMatEl> MRow;
	vector<vector<polyMatEl> > M;
	vector<vector<unsigned long> > colM(cols);								//Store the flipped version of M
	vector<polynomial> eigenPolys;										//Store the eigenvalue polynomials

	for (int i=0; i<smallMat.size(); i++)									//Create colM and M
	{
		M.push_back(MRow);
		for (int j=0; j<smallMat[i].size(); j++)
		{
			if (smallMat[i][j].isZero()==false)
			{
				M[i].push_back(polyMatEl(smallMat[i][j], j));
			}
		}
	}

	for (unsigned long i=0; i<M.size(); i++)
	{
		for (unsigned long j=0; j<M[i].size(); j++){colM[M[i][j].col].push_back(i);}			//Add entries to colM
	}

	bool allEmpty=true;											//This keeps track of when we have to stop
	for (unsigned long i=0; i<M.size(); i++)								//Check if the matrix is already empty
	{
		if (M[i].size()!=0)
		{
			allEmpty=false;
			break;
		}
	}

	while(allEmpty==false)
	{
		//printMat(M);
		//printMat(colM);
		//Look for the element with lowest (r-1)(c-1)(d+1)
		unsigned long pivR, pivC, pivCoord;
		unsigned long count=cols*M.size()*100;

		for (int i=0; i<M.size(); i++)										//Scan all elements to find the lowest count
		{
			for (int j=0; j<M[i].size(); j++)
			{
				unsigned long newCount = (M[i].size() - 1)*(colM[M[i][j].col].size() - 1)*(M[i][j].poly.degree() + 1);
				if (count > newCount)
				{
					pivR = i;
					pivC = M[i][j].col;
					count = newCount;
					pivCoord = j;
				}
			}
		}
		if (count == 0)												//Procedure for row singlets. Also functional for col singlets.
		{
			for (int i=0; i<colM[pivC].size(); i++)								//Eliminations
			{
				unsigned long curRow = colM[pivC][i];
				if (curRow!=pivR)
				{
					unsigned long index=0;
					while(M[curRow][index].col!=pivC){index++;}
					M[curRow].erase(M[curRow].begin()+index);					//Remove eliminations in M
				}
				reduce(M[curRow]);
			}
			eigenPolys.push_back(M[pivR][0].poly);								//Add to eigenPolynomials
			for (unsigned long i=0; i<M[pivR].size(); i++)							//Remove the pivot row from colM
			{
				unsigned long curCol = M[pivR][i].col;		
				unsigned long index = 0;
				while(colM[curCol][index]!=pivR){index++;}
				colM[curCol].erase(colM[curCol].begin()+index);
			}
			vector<polyMatEl>().swap(M[pivR]);								//Remove pivot from M
			vector<unsigned long>().swap(colM[pivC]);							//Remove pivot from colM
		}
		else
		{
			vector<unsigned long> pivCol = colM[pivC];							//The list of all rows we need to eliminate
			for (unsigned long i=0; i<pivCol.size(); i++)							//Start elimination
			{
				if (pivCol[i]!=pivR)
				{
					vector<polyMatEl> pivRow = M[pivR];						//Make a copy of the pivot row
					unsigned long elRow = pivCol[i];						//Name for the elimination row
					unsigned long coord=0;
					while(M[elRow][coord].col!=pivC){coord++;}					//Find sparse index of the elimination element
					polynomial pivFactor = M[elRow][coord].poly;
					polynomial elimFactor = pivRow[pivCoord].poly;
					polynomial div = gcd(pivFactor, elimFactor);					//Compute the factors we have to multiply the rows with
					pivFactor = quotient(pivFactor, div);
					pivFactor = pivFactor*(-1);
					elimFactor = quotient(elimFactor, div);

					for (unsigned long j=0; j<M[elRow].size(); j++)					//Multiply the elimination row with coefficients
					{
						M[elRow][j].poly=M[elRow][j].poly*elimFactor;
					}		

					M[elRow].reserve(M[elRow].size() + pivRow.size());				//Does this do anything?
					for (unsigned long j=0; j<pivRow.size(); j++)
					{
						unsigned long x=0;
						while (M[elRow][x].col<pivRow[j].col && x<(M[elRow].size()-1)){x++;}				//Find place for insertion
						if (M[elRow][x].col>pivRow[j].col)								//If the entry was zero
						{
							M[elRow].insert(M[elRow].begin()+x, polyMatEl(pivFactor*pivRow[j].poly, pivRow[j].col));//Insert entry
							unsigned long insertCol = pivRow[j].col;
							colM[insertCol].push_back(elRow);							//Register in the colDom M	
						}

						else if (M[elRow][x].col==pivRow[j].col)							//If the entry was nonzero
						{
							M[elRow][x].poly = M[elRow][x].poly + pivFactor*pivRow[j].poly;				//Just add to it
							if (M[elRow][x].poly.isZero())								//If this causes the entry to vanish
							{
								M[elRow].erase(M[elRow].begin()+x);
								unsigned long insertCol = pivRow[j].col;
								for (unsigned long k=0; k<colM[insertCol].size(); k++)				//Eliminate the entry from colM
								{
									if (colM[insertCol][k]==elRow)
									{
										colM[insertCol].erase(colM[insertCol].begin()+k);
										break;
									}
								}
							}
						}
			
						else 												//Special case when we need to add to the end of the row
						{
							M[elRow].push_back(polyMatEl(pivFactor*pivRow[j].poly, pivRow[j].col));			//insert at the end of the row
							unsigned long insertCol = pivRow[j].col;
							colM[insertCol].push_back(elRow);							//Register in the colDom M
						}
					}
					reduce(M[elRow]);
				}
			}

			eigenPolys.push_back(M[pivR][pivCoord].poly);							//Add the pivot to the eigenpolynomials
			for (unsigned long j=0; j<M[pivR].size(); j++)							//Remove the pivot row from colM
			{
				unsigned long curCol = M[pivR][j].col;		
				unsigned long index = 0;
				while(colM[curCol][index]!=pivR){index++;}
				colM[curCol].erase(colM[curCol].begin()+index);
			}
			vector<polyMatEl>().swap(M[pivR]);								//Remove the pivot row from M
		}
			

		allEmpty=true;
		for (unsigned long i=0; i<M.size(); i++)
		{
			if (M[i].size()!=0)
			{
				allEmpty=false;
				break;
			}
		}
	}

	vector<long> eigenvalues;										//Get all eigenvalues
	for (int i=0; i<eigenPolys.size(); i++)
	{
		if(eigenPolys[i].degree()>0)									//Solve any polynomial that's actually a polynomial 
		{
			long curContent = content(eigenPolys[i]);
			eigenPolys[i] = eigenPolys[i]/curContent;
			vector<long> solutions = solvePoly(eigenPolys[i]);
			eigenvalues.insert(eigenvalues.end(), solutions.begin(), solutions.end());
		}
	}
	return eigenvalues;
}
