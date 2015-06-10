//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The classes in this file form the computer algebra system that allows the user to input Beta-functions easily. 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "common.h"
#include "symbolics.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Variable Methods-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
bool Variable::operator <(const Variable &rhs) const
{
	if (name==rhs.name)
	{
		if (isComplex==true && rhs.isComplex==true)
		{
			if (CC==rhs.CC){return false;}
			else if (CC==true && rhs.CC==false){return false;}
			else {return true;}
		}
		else if (isComplex==false && rhs.isComplex==false){return false;}
		else if (isComplex==true && rhs.isComplex==false){return false;}
		else {return true;}
	}
	else {return (name < rhs.name);}
}	

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Param Methods----------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Param::Param(string name, unsigned long size): name(name), size(size)
{
	isComplex=false;
	CC=false;
	if (size>1)
	{
		for (unsigned long i=0; i<size; i++)
		{
			for (unsigned long j=0; j<size; j++)
			{
				string temp = name;
				ostringstream tempStringStream;
				tempStringStream << i << j;
				string tempString = tempStringStream.str(); 
				temp.append(tempString);
				vars.push_back(Variable(temp));
			}
		}
	}
	else 
	{
		vars.push_back(Variable(name));
	}
}

Param::Param(string name, unsigned long size, bool isComplex): name(name), size(size), isComplex(isComplex) 
{
	CC=false;
	if (size>1)
	{
		for (unsigned long i=0; i<size; i++)
		{
			for (unsigned long j=0; j<size; j++)
			{
				string temp = name;
				ostringstream tempStringStream;
				tempStringStream << i << j;
				string tempString = tempStringStream.str(); 
				temp.append(tempString);
				vars.push_back(Variable(temp));
			}
		}
	}

	else 
	{
		vars.push_back(Variable(name));
	}

	for (unsigned long i=0; i<vars.size(); i++){vars[i].isComplex=isComplex;}
}

	
void Param::Conjugate()
{
	if (isComplex==true)
	{
		CC = !CC;
		for (unsigned long i=0; i<vars.size(); i++){vars[i].CC=!vars[i].CC;}
	}
}

void Param::botRight()
{
	for (unsigned long i=0; i<(size*size)-1; i++){vars[i].isZero=true;}
	vars[(size*size)-1].name = name;
}

void Param::diag()
{
	for (unsigned long i=0; i<size; i++)
	{
		for (unsigned long j=0; j<size; j++)
		{
			if (i!=j){vars[i*size+j].isZero=true;}
			else 
			{
				vars[i*size+i].name = name;
				ostringstream tempStream;
				tempStream << i;
				string tempString = tempStream.str();
				vars[i*size+i].name.append(tempString);
			}
		}
	}
}
void Param::botRightDiag()
{
	for (unsigned long i=0; i<size; i++)
	{
		for (unsigned long j=0; j<size; j++)
		{
			if (i!=j){vars[i*size+j].isZero=true;}
			else 
			{
				if (i==(size-1))
				{
					vars[i*size+i].name = name;
					ostringstream tempStream;
					tempStream << size;
					string tempString = tempStream.str();
					vars[i*size+i].name.append(tempString);
				}
				else 
				{
					vars[i*size+i].name = name;	
					ostringstream tempStream;
					tempStream << 1;
					string tempString = tempStream.str();
					vars[i*size+i].name.append(tempString);			
				}
			}
		}
	}
}

void Param::print()
{
	for (unsigned long i=0; i<size; i++)
	{
		for (unsigned long j=0; j<size; j++)
		{
			if (vars[j*size+i].isZero==false)
			{
				cout << vars[j*size+i].name;
				if (vars[j*size+i].CC==true){cout << "c";}
			}
			else {cout << "0";}
			cout << " ";
				
		}
		cout << endl;		
	}
}

Function Param::operator*(const Param &B) const 
{
	if (size==1 && B.size!=1)
	//In the special case where one of the Parameters is a scalar and the other is a matrix, we treat the scalar like it is multiplied by the unity matrix. 
	//This builds that matrix and calls the method again using that matrix.
	{
		Param temp(name, B.size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			for (unsigned long j=0; j<temp.size; j++)
			{
				if (i==j){temp.vars[i*temp.size+i]=vars[0];}
				else {temp.vars[i*temp.size+j].isZero=true;}
			}
		}		
		return temp*B;
	}

	else if (size!=1 && B.size==1)
	{
		Param temp(B.name, size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			for (unsigned long j=0; j<temp.size; j++)
			{
				if (i==j){temp.vars[i*temp.size+i]=B.vars[0];}
				else {temp.vars[i*temp.size+j].isZero=true;}
			}
		}
		return (*this)*temp;
	}

	else if (size!=B.size)
	{
		cout << "non-matching error" << endl; return 0;
	}
	
	else
	{
		Function temp(size);
		for(unsigned long i=0; i<size*size; i++)
		{
			unsigned long col = i%size;
			unsigned long row = (i-col)/size;
			for (unsigned long j=0; j<size; j++)
			{
				if (vars[size*row + j].isZero==false && B.vars[B.size*j + col].isZero==false)
				{
					vector<Variable> newterm;
					newterm.push_back(vars[size*row + j]);
					newterm.push_back(B.vars[B.size*j + col]);
					temp.terms[i].push_back(newterm);
					temp.coeffs[i].push_back(ir(1));
				}
			}
		}
	return temp;
	}
}

Function Param::operator*(const Function &B) const
{
	if (B.size!=1 && size==1)
	{
		Param temp(name, B.size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			for (unsigned long j=0; j<temp.size; j++)
			{
				if (i==j){temp.vars[i*temp.size+i]=vars[0];}
				else {temp.vars[i*temp.size+j].isZero=true;}
			}
		}
		return temp*B;
	}

	else if (B.size==1 && size!=1)
	{
		Function temp(size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=B.terms[0];
			temp.coeffs[i*temp.size+i]=B.coeffs[0];
		}
		return (*this)*temp;
	}

	else if (B.size!=size)
	{
		cout << "non-matching error" << endl; return 0;
	}


	else
	{
		Function temp(B.size);
		for (unsigned long i=0; i<B.size*B.size; i++)
		{
			unsigned long col = i%B.size;
			unsigned long row = (i-col)/B.size;
			for (unsigned long j=0; j<B.size; j++)
			{
				for (unsigned long k=0; k<B.terms[row*B.size+j].size(); k++)
				{
					if (vars[size*j + col].isZero==false)
					{
						vector<Variable> newterm = B.terms[row*B.size+j][k];
						newterm.push_back(vars[size*j + col]);
						temp.terms[i].push_back(newterm);
						temp.coeffs[i].push_back(B.coeffs[row*B.size+j][k]);
					}
				}
			}
		}
		return temp;
	}
}

Function Param::operator*(const ir &r) const
{
	Function temp(size);
	for (unsigned long i=0; i<size*size; i++)
	{
		vector<Variable> newterm;
		newterm.push_back(vars[i]);
		temp.terms[i].push_back(newterm);
		temp.coeffs[i].push_back(r);
	}
	return temp;
}

Function Param::operator*(const long &a) const
{
	Function temp(size);
	for (unsigned long i=0; i<size*size; i++)
	{
		vector<Variable> newterm;
		newterm.push_back(vars[i]);
		temp.terms[i].push_back(newterm);
		temp.coeffs[i].push_back(ir(a));
	}
	return temp;
}

Param Dagger(Param A)
{
	Param temp = A;
	for (unsigned long i=0; i<temp.size;i++)
	{
		for (unsigned long j=0; j<temp.size; j++)
		{
			temp.vars[j*temp.size + i] = A.vars[i*A.size+j];
		}
	}
	temp.Conjugate();
	return temp;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Function Methods-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Function::Function(unsigned long size): size(size)
{
	vector<vector<Variable> > tempterm;
	vector<ir> temprat;
	for (unsigned long i=0; i<size*size; i++)
	{
		terms.push_back(tempterm);
		coeffs.push_back(temprat);
	}
}

void Function::print()
{
	for (unsigned long i=0; i<size; i++)
	{
		for (unsigned long j=0; j<size; j++)
		{
			cout << "(" << i << "," << j << ") = "; 
			for (unsigned long k=0; k<coeffs[i*size+j].size(); k++)
			{
				coeffs[i*size+j][k].print();
				cout << "*";
				for (unsigned long l=0; l<terms[i*size+j][k].size(); l++)
				{
					cout << terms[i*size+j][k][l].name;
					if (terms[i*size+j][k][l].CC==true){cout << "c";}
					cout << "*";
				}
				cout << " + ";
			}
			cout << endl << endl;
		}
	}
}
void Function::print(Param A)
{
	for (unsigned long i=0; i<size; i++)
	{
		for (unsigned long j=0; j<size; j++)
		{
			cout << A.vars[i*size+j].name;
			if (A.vars[i*size+j].isComplex == true && A.vars[i*size+j].CC == true){cout << "c";}
			cout <<  " = ";
			for (unsigned long k=0; k<coeffs[i*size+j].size(); k++)
			{
				coeffs[i*size+j][k].print();
				cout << "*";
				for (unsigned long l=0; l<terms[i*size+j][k].size(); l++)
				{
					cout << terms[i*size+j][k][l].name;
					if (terms[i*size+j][k][l].CC==true){cout << "c";}
					cout << "*";
				}
				cout << " + ";
			}
			cout << endl << endl;;
		}
	}
}

Function Function::operator*(const Param &B) const
{
	if (size==1 && B.size!=1)
	{
		Function temp(B.size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=terms[0];
			temp.coeffs[i*temp.size+i]=coeffs[0];
		}
		return temp*B;
	}

	else if (size!=1 && B.size==1)
	{
		Param temp(B.name, size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			for (unsigned long j=0; j<temp.size; j++)
			{
				if (i==j){temp.vars[i*temp.size+i]=B.vars[0];}
				else {temp.vars[i*temp.size+j].isZero=true;}
			}
		}
		return (*this)*temp;
	}

	else if (size!=B.size)
	{
		cout << "non-matching error" << endl; return 0;
	}


	else
	{
		Function temp(size);
		for (unsigned long i=0; i<size*size; i++)
		{
			unsigned long col = i%size;
			unsigned long row = (i-col)/size;
			for (unsigned long j=0; j<size; j++)
			{
				for (unsigned long k=0; k<terms[row*size+j].size(); k++)
				{
					if (B.vars[B.size*j + col].isZero==false)
					{
						vector<Variable> newterm = terms[row*size+j][k];
						newterm.push_back(B.vars[B.size*j + col]);
						temp.terms[i].push_back(newterm);
						temp.coeffs[i].push_back(coeffs[row*size+j][k]);
					}
				}
			}
		}
		return temp;
	}
}

Function Function::operator*(const Function &B) const
{
	if (size==1 && B.size!=1)
	{
		Function temp(B.size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=terms[0];
			temp.coeffs[i*temp.size+i]=coeffs[0];
		}
		return temp*B;
	}

	else if (size!=1 && B.size==1)
	{
		Function temp(size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=B.terms[0];
			temp.coeffs[i*temp.size+i]=B.coeffs[0];
		}
		return (*this)*temp;
	}

	else if (size!=B.size)
	{
		cout << "non-matching error" << endl; return 0;
	}

	else
	{
		Function temp(size);
		for (unsigned long i=0; i<size*size; i++)
		{
			unsigned long col = i%size;
			unsigned long row = (i-col)/size;
			for (unsigned long j=0; j<size; j++)
			{
				for (unsigned long k=0; k<terms[row*size+j].size(); k++)
				{
					for (unsigned long l=0; l<B.terms[B.size*j+col].size(); l++)
					{
						vector<Variable> newterm = terms[row*size+j][k];
						newterm.insert(newterm.end(), B.terms[B.size*j+col][l].begin(), B.terms[B.size*j+col][l].end());
						temp.terms[i].push_back(newterm);
						temp.coeffs[i].push_back(coeffs[row*size+j][k]*B.coeffs[B.size*j+col][l]);
					}
				}
			}
		}
		return temp;
	}
}

Function Function::operator*(const ir &r) const
{
	Function temp = (*this);
	for (unsigned long i=0; i<coeffs.size(); i++)
	{
		for (unsigned long j=0; j<coeffs[i].size(); j++)
		{
			temp.coeffs[i][j] = temp.coeffs[i][j]*r;
		}
	}
	return temp;
}

Function Function::operator*(const long &a) const
{
	Function temp = (*this);
	for (unsigned long i=0; i<coeffs.size(); i++)
	{
		for (unsigned long j=0; j<coeffs[i].size(); j++)
		{
			temp.coeffs[i][j] = temp.coeffs[i][j]*ir(a);
		}
	}
	return temp;
}

Function Function::operator+(const Function &B) const
{
	if (size==1 && B.size!=1)
	{
		Function temp(B.size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=terms[0];
			temp.coeffs[i*temp.size+i]=coeffs[0];
		}
		return temp+B;
	}

	else if (size!=1 && B.size==1)
	{
		Function temp(size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=B.terms[0];
			temp.coeffs[i*temp.size+i]=B.coeffs[0];
		}
		return (*this)+temp;
	}

	else if (size != B.size)
	{
		cout << "non-matching error" << endl; return 0;
	}

	else
	{
		Function temp(size);
		for (unsigned long i=0; i<size*size; i++)
		{
			temp.terms[i] = terms[i];
			temp.coeffs[i] = coeffs[i];
			for (unsigned long j=0; j<B.coeffs[i].size(); j++)
			{
				temp.terms[i].push_back(B.terms[i][j]);
				temp.coeffs[i].push_back(B.coeffs[i][j]);
			}
		}
		return temp;
	}
}

Function Function::operator-(const Function &B) const
{
	if (size==1 && B.size!=1)
	{
		Function temp(B.size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=terms[0];
			temp.coeffs[i*temp.size+i]=coeffs[0];
		}
		return temp-B;
	}

	else if (size!=1 && B.size==1)
	{
		Function temp(size);
		for (unsigned long i=0; i<temp.size; i++)
		{
			temp.terms[i*temp.size+i]=B.terms[0];
			temp.coeffs[i*temp.size+i]=B.coeffs[0];
		}
		return (*this)-temp;
	}

	else if (size != B.size)
	{
		cout << "non-matching error" << endl; return 0;
	}

	else
	{
		Function temp(size);
		for (unsigned long i=0; i<size*size; i++)
		{
			temp.terms[i] = terms[i];
			temp.coeffs[i] = coeffs[i];
			for (unsigned long j=0; j<B.coeffs[i].size(); j++)
			{
				temp.terms[i].push_back(B.terms[i][j]);
				temp.coeffs[i].push_back(ir(-1)*B.coeffs[i][j]);
			}
		}
		return temp;
	}
}

Function Tr(Function A)
{
	Function temp(1);
	for (unsigned long i=0; i<A.size; i++)
	{
		for (unsigned long j=0; j<A.coeffs[i+A.size*i].size(); j++)
		{
			temp.terms[0].push_back(A.terms[i+A.size*i][j]);
			temp.coeffs[0].push_back(A.coeffs[i+A.size*i][j]);
		}
	}
	return temp;
}

Function Tr(Param A)
{
	Function temp(1);
	for (unsigned long i=0; i<A.size; i++)
	{
		vector<Variable> tempvar;
		tempvar.push_back(A.vars[i*A.size+i]);
		temp.terms[0].push_back(tempvar);
		temp.coeffs[0].push_back(1);
	}
	return temp;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------BetaFunc Methods-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void BetaFunc::operator=(const Function &A) 
{
	func = A;
}

void BetaFunc::print()
{
	func.print(par);
}

BetaFunc Conjugate(BetaFunc A)
{
	A.par.Conjugate();
	for (unsigned long i=0; i<A.func.terms.size(); i++)
	{
		for (unsigned long j=0; j<A.func.terms[i].size(); j++)
		{
			for (unsigned long k=0; k<A.func.terms[i][j].size(); k++)
			{
				if (A.func.terms[i][j][k].isComplex==true){A.func.terms[i][j][k].CC = !A.func.terms[i][j][k].CC;}
			}
		}
	}
	return A;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------long Methods-----------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Function operator*(const long a, Param B)
{
	Function temp(B.size);
	for (unsigned long i=0; i<B.size*B.size; i++)
	{
		vector<Variable> newterm;
		newterm.push_back(B.vars[i]);
		temp.terms[i].push_back(newterm);
		temp.coeffs[i].push_back(ir(a));
	}
	return temp;
}

Function operator*(const long a, Function B)
{
	Function temp = B;
	for (unsigned long i=0; i<B.coeffs.size(); i++)
	{
		for (unsigned long j=0; j<B.coeffs[i].size(); j++)
		{
			temp.coeffs[i][j] = temp.coeffs[i][j]*ir(a);
		}
	}
	return temp;
}
