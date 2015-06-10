//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The structs in this file form the computer algebra system that allows the user to input Beta-functions easily. 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _SYMBOLICS_H_INCLUDED 
#define _SYMBOLICS_H_INCLUDED 
#include "common.h"

struct Param;												//Forward declarations
struct Function;

struct Variable
//Struct that represents a single variable.
{
	Variable(string name): name(name) {isComplex=false; CC=false; isZero=false;}
	Variable(string name, bool isComplex): name(name), isComplex(isComplex) {CC=false; isZero=false;}

	bool operator <(const Variable &rhs) const;

	string name;
	bool isComplex, CC, isZero;
};

struct Param
//Struct that represents a single parameter, which can consist of a single variable or a matrix of variables
{
	Param(string name, unsigned long size); 
	Param(string name, unsigned long size, bool isComplex);
	void Conjugate();
	void botRight();
	void diag();	
	void botRightDiag();
	void print();
	Function operator*(const Param &B) const;
	Function operator*(const Function &B) const;
	Function operator*(const ir &r) const;	
	Function operator*(const long &a) const;

	string name;
	unsigned long size;
	bool isComplex, CC;
	vector<Variable> vars;
};

Param Dagger(Param A);

struct Function
//Struct that represents a function of parameters
{
	Function(unsigned long size); 
	void print();
	void print(Param A);
	Function operator*(const Param &B) const;
	Function operator*(const Function &B) const;
	Function operator*(const ir &r) const;
	Function operator*(const long &a) const;
	Function operator+(const Function &B) const;
	Function operator-(const Function &B) const;
	
	unsigned long size; 
	vector<vector<vector<Variable> > > terms;
	vector<vector<ir> > coeffs;
};

Function Tr(Function A);
Function Tr(Param A);

struct BetaFunc
//Struct that represents a Beta-functions. Contains a Function & the corresponding parameter
{
	BetaFunc(Param par): par(par), func(Function(0)) {}
	void operator=(const Function &A);
	void print();

	Param par;
	Function func;
};

BetaFunc Conjugate(BetaFunc A);

Function operator*(const long a, Param B);
Function operator*(const long a, Function B);

#endif 
