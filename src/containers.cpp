//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//The structs in this file are the numerical containers for various objects used throughout the calculations the program performs.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
#include "common.h"
#include "containers.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Matrix Element Methods-------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
vector<MatEl> operator+(const vector<MatEl> &A, const vector<MatEl> &B)
{
	vector<MatEl> result;
	unsigned long aCoord=0, bCoord=0;
	while(aCoord!=A.size() || bCoord!=B.size())
	{
		if (aCoord==A.size())
		{
			result.push_back(B[bCoord]);
			bCoord++;
		}
		else if (bCoord==B.size())
		{
			result.push_back(A[aCoord]);
			aCoord++;
		}
		else if (A[aCoord].col<B[bCoord].col)
		{
			result.push_back(A[aCoord]);
			aCoord++;
		}
		else if (A[aCoord].col>B[bCoord].col)
		{
			result.push_back(B[bCoord]);
			bCoord++;
		}
		else if (A[aCoord].col==B[bCoord].col)
		{
			if (A[aCoord].val+B[bCoord].val!=0)
			{
				result.push_back(MatEl(A[aCoord].val+B[bCoord].val, A[aCoord].col));
			}
			aCoord++;
			bCoord++;
		}
	}
	return result;
}

void rowGCD(vector<MatEl> &R)
{
	if (R.size()==0){}
	else if (R.size()==1){R[0].val=1;}
	else
	{
		long div = gcd(R[0].val, R[1].val);
		for (unsigned long i=2; i<R.size(); i++){div = gcd(div, R[i].val);}
		for (unsigned long i=0; i<R.size(); i++){R[i].val/=div;}
	}
}

void printMat(vector<vector<MatEl> > M)
{
	for (unsigned long i=0; i<M.size(); i++)
	{
		if (M[i].size()!=0){cout << "[" << i << "] ";}
		for (unsigned long j=0; j<M[i].size(); j++)
		{
			cout << M[i][j].val << "(" << M[i][j].col << ") ";
		}
		if (M[i].size()!=0){cout << endl;}
	}
	cout << endl;
}

void printVec(vector<MatEl> V)
{
	for (unsigned long i=0; i<V.size(); i++)
	{
		cout << V[i].val << "(" << V[i].col << ") ";
	}
	cout << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Term Methods-----------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
__uint128_t termCantor(const vector<Term> &terms)
{
	deque<__uint128_t> nums;
	int MAXTERMS = 8;
	for (int i=0; i<terms.size(); i++)
	{
		nums.push_back(Cantor(makePos(terms[i].pow), terms[i].var));
	}
	for (int j=terms.size(); j<MAXTERMS; j++){nums.push_front(0);}		//Fill with zeroes until MAXTERMS

	while (true)
	{
		if (nums.size()==1){return nums.front();}
		int s = nums.size();
		for (int i=0; i<s/2; i++)
		{
			nums.push_back(Cantor(nums[0], nums[1]));
			nums.pop_front();
			nums.pop_front();
		}
	}
}

vector<Term> operator+(const vector<Term> &A, const vector<Term> &B)
{
	vector<Term> result;
	unsigned long aCoord=0, bCoord=0;
	while(aCoord!=A.size() || bCoord!=B.size())
	{
		if (aCoord==A.size())
		{
			result.push_back(B[bCoord]);
			bCoord++;
		}
		else if (bCoord==B.size())
		{
			result.push_back(A[aCoord]);
			aCoord++;
		}
		else if (A[aCoord].var<B[bCoord].var)
		{
			result.push_back(A[aCoord]);
			aCoord++;
		}
		else if (A[aCoord].var>B[bCoord].var)
		{
			result.push_back(B[bCoord]);
			bCoord++;
		}
		else if (A[aCoord].var==B[bCoord].var)
		{
			if (A[aCoord].pow+B[bCoord].pow!=0)
			{
				result.push_back(Term(A[aCoord].pow+B[bCoord].pow, A[aCoord].var));
			}
			aCoord++;
			bCoord++;
		}
	}
	return result;
}

bool compTermVec::operator()(const vector<Term> &A, const vector<Term> &B)
{
	if (A.size()!=B.size()){return A.size()<B.size();}
	else
	{
		for (unsigned long i=0; i<A.size(); i++)
		{
			if (A[i].var!=B[i].var){return A[i].var<B[i].var;}
			else if (A[i].pow!=B[i].pow){return A[i].pow<B[i].pow;}
		}
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Equation Methods-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void Equation::fillIn(vector<ir> &A)
{
	ir solu(0);
	unsigned long coef;
	for (unsigned long i=0; i<eq.size(); i++)
	{
		if (piv==eq[i].col)
		{
			coef=i;
		}
		else
		{
			solu = solu - ir(eq[i].val)*A[eq[i].col];
		}
	}
	solu = solu/(eq[coef].val);
	A[piv]=solu;
}



void printEq(vector<Equation> N)
{
	for (unsigned long i=0; i<N.size(); i++)
	{
		cout << "[" << N[i].piv << "] ";
		for (unsigned long j=0; j<N[i].eq.size(); j++)
		{
			cout << N[i].eq[j].val << "(" << N[i].eq[j].col << ") ";
		}
		cout << endl;
	}
	cout << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Container Methods------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void container::makeConstant(int dimSize)
{
	vector<vector<Term> >().swap(Powers);
	vector<long>().swap(Coefficients);

	vector<Term> empty;
	Coefficients.push_back(1);
	Powers.push_back(empty);
	for (int i=0; i<dimSize; i++){dims.push_back(0);}
}


void container::addTerm(long coeff, vector<Term> pows)
{
	Powers.push_back(pows);
	Coefficients.push_back(coeff);
}

void container::addTermSafe(long coeff, vector<Term> pows)
{
	for (int i=0; i<Powers.size(); i++)
		{
		if (pows == Powers[i])
			{
			Coefficients[i] += coeff;
			return;
			}
		}
	Powers.push_back(pows);
	Coefficients.push_back(coeff);
}

void container::print(map<int8_t, string> symbols)
{
	for (unsigned long i=0; i<monoPows.size(); i++)
	{
		if (monoPows[i].num!=0)
		{
			cout << "(" << symbols[i] << ")^";
			monoPows[i].print();
		}
	}
	cout << "[";
	for (unsigned long i=0; i<Coefficients.size(); i++)
	{
		cout << Coefficients[i];
		for (unsigned long j=0; j<Powers[i].size(); j++)
		{
			cout << "(" << symbols[Powers[i][j].var] << ")^" << static_cast<int>(Powers[i][j].pow);
		}
		cout << " + ";
	}
	cout << "dims: ";
	for (unsigned long i=0; i<dims.size(); i++) {cout << dims[i] << " ";}
	cout << endl;
}

void container::reduce()
{
	if (Coefficients.size()==1){Coefficients[0]=1;}
	else if (Coefficients.size() > 1)
	{
		long temp = gcd(Coefficients[0], Coefficients[1]);
		for (unsigned long i=2; i<Coefficients.size(); i++){temp = gcd(temp, Coefficients[i]);}
		for (unsigned long i=0; i<Coefficients.size(); i++){Coefficients[i]/=temp;}
	}
}


bool container::isMono() const
{
	if (Powers.size()==1){return true;}
	else {return false;}
}

unsigned long container::HNP() const
{
	unsigned long result = 0;
	for (int i=0; i<Powers.size(); i++)
	{
		if (Powers[i].size()>result) {result = Powers[i].size();}
	}
	return result;
}

bool container::checkDims(vector<vector<int> > dims)
{
	if (Powers.size()<2){return true;}
	else
	{
		for (int x=0; x<dims.size(); x++)
		{
			int thisDim=0;
			for (int i=0; i<Powers[0].size(); i++)
			{
				thisDim += Powers[0][i].pow * dims[x][Powers[0][i].var];
			}
			for (int i=1; i<Powers.size(); i++)
			{
				int checkDim=0;
				for (int j=0; j<Powers[i].size(); j++)
				{
					checkDim += Powers[i][j].pow * dims[x][Powers[i][j].var];
				}
				if (checkDim != thisDim) {return false;}
			}
		}
		return true;
	}
}

container container::operator*(const container &two) const
{
	unsigned long rowCount=0;
	container Out;
	for (int i=0; i<dims.size(); i++)						//Compute dims
	{
		Out.dims.push_back(dims[i]+two.dims[i]);
	}
	map<vector<Term>, long, compTermVec> termsMap;					//Compute terms
	for (unsigned long i=0; i<Coefficients.size(); i++)
	{
		for (unsigned long j=0; j<two.Coefficients.size(); j++)
		{
			long tempCoeff = Coefficients[i]*two.Coefficients[j];
			vector<Term> tempPowers = Powers[i] + two.Powers[j];

			map<vector<Term>, long>::iterator it;
			it = termsMap.find(tempPowers);

			if(it==termsMap.end())
			{
				Out.addTerm(tempCoeff, tempPowers);
				termsMap[tempPowers] = rowCount;
				rowCount++;
			}
			else
			{
				Out.Coefficients[it->second]+=tempCoeff;
			}
		}
	}
	return Out;
}

container container::pow(const long &pow) const
{
	if (pow==1){return (*this);}
	else
	{
		container Out = (*this)*(*this);
		for (long i=2; i<pow; i++){Out = Out*(*this);}
		return Out;
	}
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Monomial Methods-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void monomial::print() const
{
	cout << coeff << " ";
	if (pow!=0) {cout << "(x^" << pow << ") ";}
}

bool monomial::operator==(const monomial &m) const
{
	return pow==m.pow;
}

monomial monomial::operator*(const long &c)
{
	monomial out(coeff*c, pow);
	return out;
}

monomial monomial::operator*(const monomial &m) const			
{
	monomial out(coeff*m.coeff, pow + m.pow);
	return out;
}

monomial monomial::operator/(const monomial &m) const
{
	monomial out(coeff/m.coeff, pow-m.pow);
	return out;
}

bool monomial::operator<(const monomial &m) const
{
	return pow < m.pow;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------Polynomial Methods-----------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void polynomial::addTerm(monomial m)
{
	if (terms.size()==0){terms.push_back(m);}
	else
	{
		for (int i=0; i<terms.size(); i++)
		{
			if (m<terms[i] == false)					//Look for the first time a monomial doesnt have higher degree anymore
			{
				if (m==terms[i]) {terms[i].coeff = terms[i].coeff + m.coeff;}
				else {terms.insert(terms.begin()+i, m);}
				return;
			}
		}
		terms.push_back(m);
	}
}

void polynomial::print() const
{
	if (terms.size()==0){cout << "0";}
	else
	{
		for (int i=0; i<terms.size(); i++){terms[i].print(); cout << " + ";}
	}
}

long polynomial::degree() const
{
	if (isZero()){return 0;}
	else {return terms[0].pow;}
}

void polynomial::setZero() 
{
	vector<monomial> temp; terms=temp;
}

long polynomial::evaluate(const long &l)
{
	long solution=0;
	for (int i=0; i<terms.size(); i++)
	{
		solution += terms[i].coeff*pow(l, terms[i].pow);
	}
	return solution;
}

polynomial polynomial::operator*(const long c) const
{
	polynomial out;
	out.terms=terms;
	for (int i=0; i<out.terms.size(); i++)
	{
		out.terms[i].coeff = out.terms[i].coeff*c;
	}
	return out;
}

polynomial polynomial::operator/(const long c) const		//Assumes divisibility
{
	polynomial out;
	out.terms=terms;
	for (int i=0; i<out.terms.size(); i++)
	{
		out.terms[i].coeff = out.terms[i].coeff/c;
	}
	return out;
}

polynomial polynomial::operator*(const monomial &m) const
{
	polynomial out;
	for (int i=0; i<terms.size(); i++)
	{
		out.terms.push_back(terms[i]*m);
	}
	return out;
}

polynomial polynomial::operator*(const polynomial &p) const
{
	vector<long> coeffs(degree()+p.degree() + 1);
	for (int i=0; i<terms.size(); i++)
	{
		for (int j=0; j<p.terms.size(); j++)
		{
			coeffs[terms[i].pow + p.terms[j].pow] += terms[i].coeff*p.terms[j].coeff;
		}
	}
	polynomial out;
	for (int i=coeffs.size()-1; i>=0; i--)
	{
		if (coeffs[i]!=0) {out.terms.push_back(monomial(coeffs[i], i));}
	}
	return out;
}

polynomial polynomial::operator+(const polynomial &p) const
{
	vector<long> coeffs(max(degree(), p.degree()) + 1);
	for (int i=0; i<terms.size(); i++)
	{
		coeffs[terms[i].pow] += terms[i].coeff;
	}
	for (int i=0; i<p.terms.size(); i++)
	{
		coeffs[p.terms[i].pow] += p.terms[i].coeff;
	}
	polynomial out;
	for (int i=coeffs.size()-1; i>=0; i--)
	{
		if (coeffs[i]!=0) {out.terms.push_back(monomial(coeffs[i], i));}
	}
	return out;
}

polynomial polynomial::operator-(const polynomial &p) const							//This operator assumes the leading terms cancel
{
	vector<long> coeffs(max(degree(), p.degree()) + 1);
	for (int i=1; i<terms.size(); i++)
	{
		coeffs[terms[i].pow] += terms[i].coeff;
	}
	for (int i=1; i<p.terms.size(); i++)
	{
		coeffs[p.terms[i].pow] -= p.terms[i].coeff;
	}

	polynomial out;
	for (int i=coeffs.size()-1; i>=0; i--)
	{
		if (coeffs[i]!=0) {out.terms.push_back(monomial(coeffs[i], i));}
	}
	return out;
}

bool polynomial::operator==(const polynomial &two) const
{
	if (terms.size() != two.terms.size()){return false;}
	else
	{
		for (int i=0; i<terms.size(); i++)
		{
			if (terms[i]==two.terms[i]==false){return false;}
		}
	}
	return true;
}


void print (const vector<polynomial> &row)
{
	for (int i=0; i<row.size(); i++)
	{
		if (row[i].isZero()==false)
		{
			cout << i << ": ";
			row[i].print();
			cout << " | ";
		}
	}
	cout << endl;
}

bool isEqual(const vector<polynomial> &A, const vector<polynomial> &B)
{
	for (int i=0; i<A.size(); i++)
	{
		if ((A[i]==B[i])==false){return false;}
	}
	return true;
}

void reduceNum(vector<polynomial> &row)   
{
	queue<long> nums;
	for (int i=0; i<row.size(); i++)
	{
		for (int j=0; j<row[i].terms.size(); j++)
		{
			nums.push(row[i].terms[j].coeff);
		}
	}

	long div = 1;
	if (nums.size()==1) {div = nums.front();}
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
		for (int j=0; j<row[i].terms.size(); j++)
		{
			row[i].terms[j].coeff = row[i].terms[j].coeff/div;
		}
	}
}

//----------------------------------------------------Functions that are used in polynomial gcd computations----------------------------------
long content(const polynomial &f)
{
	queue<long> nums;
	for (int i=0; i<f.terms.size(); i++)
	{
		nums.push(f.terms[i].coeff);
	}

	long div = 1;
	if (nums.size()==0){return 1;}
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
	return div;
}


polynomial prem(const polynomial &f, const polynomial &g)
{
	polynomial r = f*pow(g.LC(), f.degree() - g.degree() + 1);
	while(r.isZero()==false && r.degree() >= g.degree())
	{
		monomial u = r.LT()/g.LT();
		r = r - g*u;
	}
	return r;
}

polynomial quotient(const polynomial &f, const polynomial &g)
{
	polynomial q;
	polynomial r = f;
	while(r.isZero()==false && r.degree() >= g.degree())
	{
		monomial u = r.LT()/g.LT();
		q.terms.push_back(u);
		r = r - g*u;
	}
	return q;
}


polynomial gcd(const polynomial &f, const polynomial &g)
{
	vector<polynomial> prems;
	if (f.degree()<g.degree())
	{
		prems.push_back(g);
		prems.push_back(f);
	}
	else
	{
		prems.push_back(f);
		prems.push_back(g);
	}
	while(prems.back().isZero()==false)
	{
		prems.push_back(prem(prems[prems.size()-2], prems[prems.size()-1]));
		long curContent = content(prems.back());
		prems.back() = prems.back()/curContent;
	}
	polynomial result = prems[prems.size()-2];
	if (prems.size()==3)										//If g is the gcd, make sure it is primitivised
	{
		long curContent = content(result);
		result = result/curContent;
	}
	long fgcontent = gcd(content(f), content(g));
	result = result*fgcontent;
	return result;
}

