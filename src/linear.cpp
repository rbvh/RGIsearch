//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This file contains the functions that handle the solving of the linear systems in regular polynomial searching
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "linear.h"
#include "common.h"
#include "containers.h"

void elimOneRow(vector<vector<MatEl> > &M, vector<vector<unsigned long> > &colM, queue<unsigned long> &elimCandidates, vector<Equation> &N)
//Elimination functions for singlet rows/cols
{
	while(elimCandidates.empty()==false)
	{
		unsigned long elimRow = elimCandidates.front();
		elimCandidates.pop();
							
		if (M[elimRow].size()!=0)
		{
			unsigned long elimCol = M[elimRow][0].col;
			N.push_back(Equation(M[elimRow], elimCol));
			vector<MatEl>().swap(M[elimRow]);							//Remove elimination row 
	
			vector<unsigned long> removeRows = colM[elimCol];
			vector<unsigned long>().swap(colM[elimCol]);						//Remove elimination col
		
			for (unsigned long i=0; i<removeRows.size(); i++)
			{
				unsigned long curRow = removeRows[i];
				if (curRow!=elimRow)
				{
					unsigned long j=0;
					while(M[curRow][j].col!=elimCol){j++;}
					M[curRow].erase(M[curRow].begin()+j);
					if (M[curRow].size()==1)
					{
						elimCandidates.push(curRow);
					}
				}
			}
		}
	}
}
		
void elimMoreRow(vector<vector<MatEl> > &M, vector<MatEl> pivRow, vector<vector<unsigned long> > &rowC, vector<vector<unsigned long> > &colC, vector<vector<unsigned long> > &colM, unsigned long elR, unsigned long pivC, unsigned long pivCoord)
//Elimination function for multiplet rows/cols
{
	unsigned long coord=0;
	while(M[elR][coord].col!=pivC){coord++;}								//Find sparse index of the elimlongion element
	long div = gcd(M[elR][coord].val, pivRow[pivCoord].val);						//Fraction free elimlongion
	long pivCoef = -M[elR][coord].val/div;
	long elCoef = pivRow[pivCoord].val/div;
	for (unsigned long j=0; j<M[elR].size(); j++){M[elR][j].val*=elCoef;}					//Multiply row with coefficients
	unsigned long insertRowCount = M[elR].size();								//Store size of the row
	M[elR].reserve(M[elR].size() + pivRow.size());
	for (unsigned long j=0; j<pivRow.size(); j++)
	{
		unsigned long x=0;
		while (M[elR][x].col<pivRow[j].col && x<(M[elR].size()-1)){x++;}				//Find place for insertion
		if (M[elR][x].col>pivRow[j].col)								//If the entry was zero
		{
			M[elR].insert(M[elR].begin()+x, MatEl(pivCoef*pivRow[j].val, pivRow[j].col));		//Insert entry
			unsigned long insertCol = pivRow[j].col;
			unsigned long insertColCount = colM[insertCol].size();
			colM[insertCol].push_back(elR);								//Register in the colDom M
			for (unsigned long k=0; k<colC[insertColCount].size(); k++)				//+1 to the colCount of the current col
			{
				if (colC[insertColCount][k]==insertCol)
				{
					colC[insertColCount].erase(colC[insertColCount].begin()+k);
					colC[insertColCount].shrink_to_fit();
					break;
				}
			}
			colC[insertColCount+1].push_back(insertCol);	
		}

		else if (M[elR][x].col==pivRow[j].col)								//If the entry was nonzero
		{
			M[elR][x].val += pivCoef*pivRow[j].val;							//Just add to it
			if (M[elR][x].val==0)									//If this causes the entry to vanish
			{
				M[elR].erase(M[elR].begin()+x);
				unsigned long insertCol = pivRow[j].col;
				unsigned long insertColCount = colM[insertCol].size();
				for (unsigned long k=0; k<colM[insertCol].size(); k++)				//Eliminate the entry from colM
				{
					if (colM[insertCol][k]==elR)
					{
						colM[insertCol].erase(colM[insertCol].begin()+k);
						colM[insertCol].shrink_to_fit();
						break;
					}
				}

				for (unsigned long k=0; k<colC[insertColCount].size(); k++)			//-1 to the colCount of the current col
				{
					if (colC[insertColCount][k]==insertCol)
					{
						colC[insertColCount].erase(colC[insertColCount].begin()+k);
						colC[insertColCount].shrink_to_fit();
						break;
					}
				}
				colC[insertColCount-1].push_back(insertCol);
			}
		}
			
		else 												//Special case when we need to add to the end of the row
		{
			M[elR].push_back(MatEl(pivCoef*pivRow[j].val, pivRow[j].col));				//insert at the end of the row
			unsigned long insertCol = pivRow[j].col;
			unsigned long insertColCount = colM[insertCol].size();
			colM[insertCol].push_back(elR);								//Register in the colDom M
			for (unsigned long k=0; k<colC[insertColCount].size(); k++)				//+1 to the colCount of the current col
			{
				if (colC[insertColCount][k]==insertCol)
				{
					colC[insertColCount].erase(colC[insertColCount].begin()+k);
					colC[insertColCount].shrink_to_fit();
					break;
				}
			}
			colC[insertColCount+1].push_back(insertCol);
		}
	}

	rowGCD(M[elR]);
	M[elR].shrink_to_fit();

	for(unsigned long k=0; k<rowC[insertRowCount].size(); k++)						//Remove RowCount of the elimination row
	{
		if (rowC[insertRowCount][k]==elR)
		{
			rowC[insertRowCount].erase(rowC[insertRowCount].begin()+k);
			rowC[insertRowCount].shrink_to_fit();
			break;
		}
	}
	rowC[M[elR].size()].push_back(elR);									//And insert the new one

}
//------------------------------------------------------------------------------------------------------------------------------------------------
void OneCounts(vector<vector<MatEl> > &M, unsigned long cols, vector<Equation> &N)
//Eliminate all row singlets
{
	queue<unsigned long> elimCandidates;
	vector<vector<unsigned long> > colM(cols);
	for (unsigned long i=0; i<M.size(); i++)
	{
		for (unsigned long j=0; j<M[i].size(); j++){colM[M[i][j].col].push_back(i);}			//Add entries to colM
		if (M[i].size()==1){elimCandidates.push(i);}							//Find initial possible rowcount candidates
	}

	elimOneRow(M, colM, elimCandidates, N);

	vector<vector<MatEl> > newM;
	for (unsigned long i=0; i<M.size(); i++)
	{
		if (M[i].size()!=0){newM.push_back(M[i]);}
	}
	M=newM;
}
//----------------------------------------------------------------Performs Markowitz Pivoting-------------------------------------------------------
void Markowitz(vector<vector<MatEl> > &M, unsigned long cols, vector<Equation> &N)
{
//Perform Gaussian Elimination with Markowitz Pivoting
	OneCounts(M, cols, N);											//First eliminate all row singlets
	vector<vector<unsigned long> > colM(cols);
	vector<vector<unsigned long> > rowC(cols+1);
	vector<vector<unsigned long> > colC(M.size()+1);
	for (unsigned long i=0; i<M.size(); i++)
	{
		for (unsigned long j=0; j<M[i].size(); j++){colM[M[i][j].col].push_back(i);}			//Add entries to colM
		if (M[i].size()!=0) {rowC[M[i].size()].push_back(i);}
	}

	for (unsigned long i=0; i<colM.size(); i++)								//Create colC
	{
		if (colM[i].size()!=0) {colC[colM[i].size()].push_back(i);}
	}

	bool allEmpty=true;
	for (unsigned long i=0; i<M.size(); i++)
	{
		if (M[i].size()!=0)
		{
			allEmpty=false;
			break;
		}
	}

	unsigned long m=0;
	while(allEmpty==false)											//Loop over eliminations
	{
		unsigned long pivR, pivC, pivCoord=0; 
		unsigned long count = 2000000000;								//Set count to large number

		if (rowC[1].size()!=0 || colC[1].size()!=0)							//Checking for singleton rows/cols
		{
			if (rowC[1].size()!=0)
			{
				pivR = rowC[1][0];
				pivC = M[pivR][0].col;
				count = 0;
			}
	
			else 
			{
				pivC = colC[1][0];
				pivR = colM[pivC][0];
				count = 0;
			}
		}

		else 												//Minimize count
		{
			unsigned long i=2;
			while(true)
			{
				if ((i-1)*(i-1)>count || i>=rowC.size() || i>=colC.size()){break;}		//Stopping critereum
				for (unsigned long j=0; j<rowC[i].size(); j++)	//Search the rows
				{
					unsigned long searchR = rowC[i][j];
					for (unsigned long k=0; k<i; k++)
					{
						if ((i-1)*(colM[M[searchR][k].col].size()-1) < count)
						{
							pivR = searchR;
							pivC = M[searchR][k].col;
							count = (i-1)*(colM[M[searchR][k].col].size()-1);
						}	
					}
				}

				for (unsigned long j=0; j<colC[i].size(); j++)	//Search the cols
				{
					unsigned long searchC = colC[i][j];
					for (unsigned long k=0; k<i; k++)
					{
						if ((i-1)*(M[colM[searchC][k]].size()-1) < count)
						{
							pivC = searchC;
							pivR = colM[searchC][k];
							count = (i-1)*(M[colM[searchC][k]].size()-1);
						}
					}
				}
				i++;
			}
		}
		while(M[pivR][pivCoord].col!=pivC){pivCoord++;}							//Find the sparse index of the pivot element
		//cout << "The Count is " <<  count << " at row " << pivR << " and col " << pivC << endl;

		vector<unsigned long> pivCol = colM[pivC];
		vector<MatEl> pivRow = M[pivR];
		for (unsigned long i=0; i<pivCol.size(); i++)							//Start elimination
		{
			if (pivCol[i]!=pivR)
			{
				elimMoreRow(M, pivRow, rowC, colC, colM, pivCol[i], pivC, pivCoord);
			}
		}
		
		for (unsigned long j=0; j<M[pivR].size(); j++)							//Remove the pivot row from all structures
		{
			unsigned long curCol = M[pivR][j].col;
			unsigned long curColCount = colM[curCol].size();
			for (unsigned long k=0; k<colC[curColCount].size(); k++)				//Remove all column count entries
			{
				if (colC[curColCount][k]==curCol)
				{
					colC[curColCount].erase(colC[curColCount].begin()+k);
					colC[curColCount].shrink_to_fit();
					break;
				}
			}

			if (curColCount!=1)									//Add lowered col counts if we werent in the eliminated col
			{
				colC[curColCount-1].push_back(curCol);
			}

			for (unsigned long k=0; k<colM[curCol].size(); k++)					//Modify colDom M
			{
				if (colM[curCol][k]==pivR)
				{
					colM[curCol].erase(colM[curCol].begin()+k);
					colM[curCol].shrink_to_fit();
					break;
				}
			}
		}
		unsigned long curRowCount = M[pivR].size();
		for (unsigned long k=0; k<rowC[curRowCount].size(); k++)
		{
			if (rowC[curRowCount][k]==pivR)
			{
				rowC[curRowCount].erase(rowC[curRowCount].begin()+k);
				rowC[curRowCount].shrink_to_fit();
				break;
			}
		}									
		N.push_back(Equation(M[pivR], pivC));
		vector<MatEl>().swap(M[pivR]);

		m++;
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
} 

//------------------------------------------------------------------------------------------------------------------
vector<vector<long> > findNull(vector<vector<MatEl> > &M, unsigned long cols)
//Compute the nullspace of M using the above functions
{
	vector<Equation> N;
	vector<bool> detVars;
	for (unsigned long i=0; i<cols; i++){detVars.push_back(false);}

	Markowitz(M,cols,N);

	for (unsigned long i=0; i<N.size(); i++){detVars[N[i].piv]=true;}
	vector<vector<ir> > nullSpace;
	vector<ir> temp;
	for (unsigned long i=0; i<cols; i++){temp.push_back(ir(0));}
	for (unsigned long i=0; i<cols; i++)
	{
		if (detVars[i]==false)
		{
			nullSpace.push_back(temp);
			nullSpace[nullSpace.size()-1][i]=ir(1);
		}
	}

	for (unsigned long i=0; i<nullSpace.size(); i++)
	{
		if (N.size()!=0)
		{
			for (long j=N.size()-1; j>=0; j--)			
			{
				N[j].fillIn(nullSpace[i]);
			}
		}
	}

	vector<vector<long> > longNullSpace;
	vector<long> longTemp;
	for (unsigned long i=0; i<nullSpace.size(); i++)
	{
		longNullSpace.push_back(longTemp);
		long totalDenum=nullSpace[i][0].den;
		for (unsigned long j=1; j<nullSpace[i].size(); j++)
		{
			if (totalDenum%(nullSpace[i][j].den)!=0){totalDenum*=nullSpace[i][j].den;}
		}
		for (unsigned long j=0; j<nullSpace[i].size(); j++)
		{
			nullSpace[i][j] = nullSpace[i][j]*ir(totalDenum);
			longNullSpace[i].push_back(nullSpace[i][j].num);
		}
	}
	return longNullSpace;
}

//---------------------------------------------------------------------------------------------------------------------
unsigned long findNulldim(vector<vector<MatEl> > &M, unsigned long cols)
//Compute just the dimension of the nullspace using the above functions
{
	vector<Equation> N;
	Markowitz(M,cols,N);
	return (cols-N.size());
}

