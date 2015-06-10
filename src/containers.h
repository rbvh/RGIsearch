//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The structs in this file are the numerical containers for various objects used throughout the calculations the program performs.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef _CONTAINERS_H_INCLUDED 
#define _CONTAINERS_H_INCLUDED
#include "common.h"

struct MatEl
{
//Struct that contains sparse matrix elements. Matrices are stored in row-dominant form, meaning the row structure is 
//maintained and the nonzero elements are characterized by their value and their column index.
	MatEl(long val, unsigned long col): val(val), col(col) {}
	long val;
	unsigned long col;
};

	vector<MatEl> operator+(const vector<MatEl> &A, const vector<MatEl> &B);
	void rowGCD(vector<MatEl> &R);
	void printMat(vector<vector<MatEl> > M);
	void printVec(vector<MatEl> V);

struct Term
{
//Struct that represents a variable to some power in the computational part of the program. It has an identifier for the variable and the power of that variable
	Term(int pow, int var):var(var), pow(pow) {}
	int8_t var;
	int16_t pow;

	bool operator<(const Term &t) const {return var < t.var;}
	bool operator==(const Term &t) const {return (var==t.var && pow==t.pow);}
};

__uint128_t termCantor(const vector<Term> &terms);
vector<Term> operator+(const vector<Term> &A, const vector<Term> &B);

struct compTermVec
{
//Struct that handles comparison between vectors of terms. Used for mapping purposes
	bool operator()(const vector<Term> &A, const vector<Term> &B);
};

struct Equation
//Struct that stores equations that are eliminated in the Gaussian Elimination procedure
{
	Equation(vector<MatEl> eq, unsigned long piv): eq(eq), piv(piv) {}
	void fillIn(vector<ir> &A);

	vector<MatEl> eq;
	unsigned long piv;
};

struct container
{
//Container that stores the numerical form of Beta functions and invariants. Used throughout the program.
	container() {}
	container(const vector<int> &newDims) {dims = newDims;}

	void makeConstant(int dimSize);
	void addTerm(long coeff, vector<Term> pows);
	void addTermSafe(long coeff, vector<Term> pows);
	void print(map<int8_t, string> symbols);	
	void reduce();	
	bool isMono() const;
	unsigned long HNP() const;	
	bool checkDims(vector<vector<int> > dims);
	container operator*(const container &two) const;
	container pow(const long &pow) const;
	
	vector<int> dims;
	vector<long> Coefficients;
	vector<vector<Term> > Powers;
	vector<ir> monoPows;
};

struct monomial
//Struct that holds a monomial in the procedure of factorized searching. 
//Since the program currently only supports the factorization of a single variable, these monomials only contain a single variable, making them virutally identical to the Term class.
//However, they are separated for future improvement.
{
	monomial(long c, long p): coeff(c), pow(p) {}
	monomial(): coeff(1) {}

	void print() const;
	void setPow(long p){pow = p;}
	void setCoeff(long c){coeff = c;}
	bool operator==(const monomial &m) const;
	monomial operator*(const long &c);
	monomial operator*(const monomial &m) const;			
	monomial operator/(const monomial &m) const;
	bool operator<(const monomial &m) const;

	long coeff;
	long pow;
};

struct polynomial
//Struct that holds multiple monomial terms, e.g. a polynomial. Used in factorized searching.
{
	polynomial(){}
	polynomial(const monomial &firstTerm){terms.push_back(firstTerm);}

	void addTerm(monomial m);
	void print() const;
	monomial LT() {return terms[0];}							//Return the variable of the leading term of the polynomial
	monomial LT() const {return terms[0];} 
	long LC() const {return terms[0].coeff;}						//Return the coefficient of the leading term
	monomial LM() const {monomial out = terms.front(); out.setCoeff(1); return out;}	//Return entire leading monomial
	bool isZero() const {return terms.size()==0;}
	long degree() const;
	void setZero();
	long evaluate(const long &l);

	polynomial operator*(const long c) const;
	polynomial operator/(const long c) const;						//Assumes divisibility
	polynomial operator*(const monomial &m) const;
	polynomial operator*(const polynomial &p) const;
	polynomial operator+(const polynomial &p) const;
	polynomial operator-(const polynomial &p) const;					//This operator assumes the leading terms cancel
	bool operator==(const polynomial &two) const;

	vector<monomial> terms;
};

void print (const vector<polynomial> &row);
bool isEqual(const vector<polynomial> &A, const vector<polynomial> &B);
void reduceNum(vector<polynomial> &row);  
long content(const polynomial &f);
polynomial prem(const polynomial &f, const polynomial &g);
polynomial quotient(const polynomial &f, const polynomial &g);
polynomial gcd(const polynomial &f, const polynomial &g);

struct polyMatEl
//Sparse matrix element that contains a polynomial instead of a number. Used in factorized searching.
{
	polyMatEl(polynomial poly, unsigned long col): poly(poly), col(col) {}
	unsigned long col;
	polynomial poly;
};
#endif
