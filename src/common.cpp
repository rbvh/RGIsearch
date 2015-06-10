#include "common.h"
#include "symbolics.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------gcd Function-----------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
long gcd(long a, long b)
//Compute the gcd of two integers
{
	if (a==0 && b==0) {return 1;}
	long temp;
	while (b!=0)
	{
		temp = b;
		b = a%b;
		a=temp;
	}
	return a;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------ir Methods-------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
ir ir::operator+(const ir &r) const
{
	long tempn = this->num*r.den + r.num*this->den;
	long tempd = this->den*r.den;
	long g = gcd(tempn, tempd);
	return ir(tempn/g, tempd/g);  
}

ir ir::operator-(const ir &r) const 
{
	long tempn = this->num*r.den - r.num*this->den;
	long tempd = this->den*r.den;
	long g = gcd(tempn, tempd);
	return ir(tempn/g, tempd/g); 
}


ir ir::operator*(const ir &r) const 
{
	long tempn = this->num*r.num;
	long tempd = this->den*r.den;
	long g = gcd(tempn, tempd);
	return ir(tempn/g, tempd/g); 
}

ir ir::operator/(const ir &r) const
{
	long tempn = this->num*r.den;
	long tempd = this->den*r.num;
	long g = gcd(tempn, tempd);
	return ir(tempn/g, tempd/g);
}


ir ir::IRgcd(const ir &b) const
{
	return ir(gcd(num*b.den, den*b.num), den*b.den);
}


void ir::operator +=(const ir &r)
{
	long tempn = this->num*r.den + r.num*this->den;
	long tempd = this->den*r.den;
	long g = gcd(tempn, tempd);
	num = tempn/g;
	den = tempd/g;
}

void ir::operator -=(const ir &r)
{
	long tempn = this->num*r.den - r.num*this->den;
	long tempd = this->den*r.den;
	long g = gcd(tempn, tempd);
	num = tempn/g;
	den = tempd/g;
}

bool ir::operator!=(const ir &r) const 
{
	return (num!=r.num) || (den!=r.den);
}

void ir::print() const
{
	if (den==1){cout << num;}
	else {cout << "(" << num << "/" << den << ")";}
}

Function ir::operator*(const Param &B) const
{
	Function temp(B.size);
	for (unsigned long i=0; i<B.size*B.size; i++)
	{
		vector<Variable> newterm;
		newterm.push_back(B.vars[i]);
		temp.terms[i].push_back(newterm);
		temp.coeffs[i].push_back((*this));
	}
	return temp;
}

Function ir::operator*(const Function &B) const
{
	Function temp = B;
	for (unsigned long i=0; i<B.coeffs.size(); i++)
	{
		for (unsigned long j=0; j<B.coeffs[i].size(); j++)
		{
			temp.coeffs[i][j] = temp.coeffs[i][j]*(*this);
		}
	}
	return temp;
}

ir gcd(ir a, ir b)
{
	return a.IRgcd(b);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Cantor Pairing Functions-----------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------


__uint128_t Cantor(const __uint128_t &a, const __uint128_t &b)
//Represents the Cantor pairing function: A bijection from two integers to one. Used to produce unique integers for a monomial term to form a more efficient database 
{
	return (a + b) * (a + b + 1) / 2 + a;
}

long Cantor(const long &a, const long &b) 
{
	return (a + b) * (a + b + 1) / 2 + a;
}

__uint128_t makePos(const int &a)
//Helper function for the procedure of using the Cantor pairing function for the database
{
	return a>=0 ? 2*a : -2*a-1;
}
