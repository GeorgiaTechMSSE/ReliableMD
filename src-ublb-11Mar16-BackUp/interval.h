// geninterval.h: interface for the interval class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GENERALIZED_INTERVAL_H_)
#define _GENERALIZED_INTERVAL_H_


#include "IVException.h"
// #include <stdio.h>
#include <iostream>
using namespace std;

//enum BOOLEAN { FALSE, TRUE };

struct intervalItem		// store lower bound, upper bound values only
{
   double inf;
   double sup;
};


class interval  
{
	friend class PairEAMRad;
	friend class Pair;
	friend class ComputePressure;
	friend class ComputeTemp;
private:
	double inf;
	double sup;

public:
	static const double POS_INF;
	static const double NEG_INF;

public:
	//Arithmetic operations
	interval operator+(const interval& right);
	interval operator-(const interval& right);
	interval operator*(const interval& right);
	interval operator/(const interval& right);
	interval operator+(const double right);
	interval operator-(const double right);
	interval operator*(const double right);
	interval operator/(const double right);
	interval& operator=(const interval& right);
	interval& operator+=(const interval& right);
	interval& operator-=(const interval& right);
	interval& operator*=(const interval& right);
	interval& operator/=(const interval& right);
	interval& operator=(const double right);
	interval& operator+=(const double right);
	interval& operator-=(const double right);
	interval& operator*=(const double right);
	interval& operator/=(const double right);
	

    //		  operator int() const;
//    	operator double() const;
    
	interval dual();
	interval prop();
	interval impr();
    intervalItem item();

	//mathmatical operations/functions
	friend interval operator+(double left, const interval& right);
	friend interval operator-(double left, const interval& right);
	friend interval operator*(double left, const interval& right);
	friend interval operator/(double left, const interval& right);
	friend interval operator-(const interval& right);
	friend interval sin(const interval& p);
	friend interval cos(const interval& p);
	
	// Operation that requires other fundamental ops to be defined
	/*interval power(double exp);
    interval sqrt(interval v);*/

	//Set operations
	interval union_with(const interval& right);
	interval intersect_with(const interval& right);

	//relations
	int operator==(const interval &right);
	int operator>=(const interval &right);
	int operator>(const interval &right);
	int operator<=(const interval &right);
	int operator<(const interval &right);
	int operator==(const double right);
	int operator>=(const double right);
	int operator>(const double right);
	int operator<=(const double right);
	int operator<(const double right);

	int EQ(const interval &right);
	int GT_EQ(const interval &right);
	int GT(const interval &right);
	int LT_EQ(const interval &right);
	int LT(const interval &right);
	int S_GT_EQ(const interval &right);
	int S_GT(const interval &right);
	int S_LT_EQ(const interval &right);
	int S_LT(const interval &right);
	int INCLUDE(const interval &right);

//    size_t sizeof();

	//Modifiers
//	void setNominal(double nominal);
	void setInf(double L);
	void setSup(double R);
	void setLBound(double l);
	void setUBound(double u);

	//Access methods
	double Inf();
	double Sup();
	double LBound();
	double UBound();
//	double Nominal();
	double Mid();
	double Width();
	int			isProperinterval();
	int			isImproperinterval();
	int			isPointinterval();
	interval	SET();

	//copy constructor
	interval(const interval &src);

	//constructors and destructors
	interval(double value);
	interval(intervalItem value);
	interval(double l, double u);
	interval(double l, double n, double u);
	interval();
	virtual ~interval();

	friend ostream &operator<<(ostream &stream, const interval &itv);
};

// print_interval()
// Print an "interval object"

/*void interval::print_interval(char *hdr)
{
	printf(" interval %s is [ %f, %f] \n", hdr,inf,sup);
}*/

#ifndef _CONSTANT_PI_
#define _CONSTANT_PI_
const double PI = 3.1415926535897932384626;
#endif


#endif // !defined(_GENERALIZED_INTERVAL_H_)
/* -------------------------------------------- */