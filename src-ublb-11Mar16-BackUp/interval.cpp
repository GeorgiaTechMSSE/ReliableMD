// interval.cpp: implementation of the interval class.
//
#include "interval.h"
// #include <stdio.h>
#include <math.h>

const double interval::POS_INF = 1e200;
const double interval::NEG_INF = -1e200;


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

interval::interval()
{
	inf = 0.0;
//	nominal = 0.0;
	sup = 0.0;
}

interval::interval(double l, double r)
{
//	if(l>u)	throw IVException("IVException: Construct empty interval.");
//	nominal = (l+r)/2;
	inf = l;
	sup = r;
}
/*
interval::interval(double l, double n, double r)
{
	if(n<l && n<r || n>l && n>r)
	{
		cerr<<"["<<l<<", "<<n<<", "<<r<<"]"<<endl;
		throw IVException("IVException: Construct invalid interval - nominal value.");
	}
	inf = l;
	nominal = n;
	sup = r;
}
*/
interval::interval(double value)
{
	inf = value;
//	nominal = value;
	sup = value;
}

interval::interval(intervalItem value)
{
	inf = value.inf;
	sup = value.sup;
}

interval::~interval()
{

}

/////////////////////
// copy constructor
/////////////////////
interval::interval(const interval &src)
{
	inf = src.inf;
//	nominal = src.nominal;
	sup = src.sup;
}

//////////////////////////////////////////////////////////////////////
// Access functions
//////////////////////////////////////////////////////////////////////

double interval::Sup()
{
	return sup;
}

double interval::Inf()
{
	return inf;
}
/*
double interval::Nominal()
{
	return nominal;
}
*/
double interval::Mid()
{
	return (sup+inf)/2.0;
}

double interval::Width()
{
	return fabs(sup-inf);
}

int interval::isProperinterval()
{
	return inf<=sup ? 1 : 0;
}

int interval::isImproperinterval()
{
	return inf>=sup ? 1 : 0;
}

int interval::isPointinterval()
{
	return inf==sup ? 1 : 0;
}

interval interval::SET()
{
//	return interval(this->LBound(), this->Nominal(), this->UBound());
	return interval(this->LBound(), this->UBound());
}

double interval::UBound()
{
	return inf >= sup ? inf : sup;
}

double interval::LBound()
{
	return inf <= sup ? inf : sup;
}

//////////////////////////////////////////////////////////////////////
// Modifiers
//////////////////////////////////////////////////////////////////////

void interval::setSup(double r)
{
	sup = r;
}

void interval::setInf(double l)
{
	inf = l;
}

void interval::setUBound(double r)
{
	sup = r;
}

void interval::setLBound(double l)
{
	inf = l;
}
/*
void interval::setNominal(double nominal)
{
	this->nominal = nominal;
}
*/

//////////////////////////////////////////////////////////////////////
// I/O overload
//////////////////////////////////////////////////////////////////////

ostream &operator <<(ostream &stream, const interval &itv)
{
//	stream << "[" << itv.inf <<", "<< itv.nominal <<", "<< itv.sup << "]";
	stream << "[" << itv.inf <<", "<< itv.sup << "]";
	return stream;
}


//////////////////////////////////////////////////////////////////////
// Operator overloads
//////////////////////////////////////////////////////////////////////

interval interval::operator +(const interval &right)
{
	interval result = right;
//	if(result.isEmpty())	
//		throw IVException("IVException: add an empty interval.");
	result.setInf(inf+right.inf);
//	result.setNominal(nominal+right.nominal);
	result.setSup(sup+right.sup);
	return result;
}

interval interval::operator -(const interval &right)
{
	interval result = right;
//	if(result.isEmpty())	
//		throw IVException("IVException: subtract an empty interval.");
	result.setInf(inf-right.sup);
//	result.setNominal(nominal-right.nominal);
	result.setSup(sup-right.inf);
	return result;
}

interval interval::operator *(const interval &right)
{
//	interval result = right;
//	if(result.isEmpty())	
//		throw IVException("IVException: multiply an empty interval.");
	double a1b1, a1b2, a2b1, a2b2;
	a1b1 = inf*right.inf;
	a1b2 = inf*right.sup;
	a2b1 = sup*right.inf;
	a2b2 = sup*right.sup;
/*	int min = 0, max = 0;
	for(int i=1; i<4; i++)
	{
		if(m[i]<m[min]) min=i;
		if(m[i]>m[max]) max=i;
	}
	result.setL(m[min]);
	result.setNominal(nominal*right.nominal);
	result.setR(m[max]);
	return result;
*/
	if ( inf>=0 && sup>=0 && right.inf>=0 && right.sup>=0 )
		return interval(a1b1, a2b2);
	else if ( inf>=0 && sup>=0 && right.inf>=0 && right.sup<0 )
		return interval(a1b1, a1b2);
	else if ( inf>=0 && sup>=0 && right.inf<0 && right.sup>=0 )
		return interval(a2b1, a2b2);
	else if ( inf>=0 && sup>=0 && right.inf<0 && right.sup<0 )
		return interval(a2b1, a1b2);
	else if ( inf>=0 && sup<0 && right.inf>=0 && right.sup>=0 )
		return interval(a1b1, a2b1);
	else if ( inf>=0 && sup<0 && right.inf>=0 && right.sup<0 )
		return interval((a2b2>=a1b1?a2b2:a1b1), (a2b1<=a1b2?a2b1:a1b2));
	else if ( inf>=0 && sup<0 && right.inf<0 && right.sup>=0 )
		return interval(0.0,0.0);
	else if ( inf>=0 && sup<0 && right.inf<0 && right.sup<0 )
		return interval(a2b2, a1b2);
	else if ( inf<0 && sup>=0 && right.inf>=0 && right.sup>=0 )
		return interval(a1b2, a2b2);
	else if ( inf<0 && sup>=0 && right.inf>=0 && right.sup<0 )
		return interval(0.0, 0.0);
	else if ( inf<0 && sup>=0 && right.inf<0 && right.sup>=0 )
		return interval((a1b2<=a2b1?a1b2:a2b1), (a1b1>=a2b2?a1b1:a2b2));
	else if ( inf<0 && sup>=0 && right.inf<0 && right.sup<0 )
		return interval(a2b1, a1b1);
	else if ( inf<0 && sup<0 && right.inf>=0 && right.sup>=0 )
		return interval(a1b2, a2b1);
	else if ( inf<0 && sup<0 && right.inf>=0 && right.sup<0 )
		return interval(a2b2, a2b1);
	else if ( inf<0 && sup<0 && right.inf<0 && right.sup>=0 )
		return interval(a1b2, a1b1);
	else
		return interval(a2b2, a1b1);
	
}

interval interval::operator/(const interval& right)
{
	interval a = right;
	double a1b1, a1b2, a2b1, a2b2;

//	if(a.isEmpty())	
//		throw IVException("IVException: divide an empty interval.");
	// dividend does not include zero
/*	if(!a.INCLUDE(0))
	{
		interval div(1.0/right.R, 1.0/right.nominal, 1.0/right.inf);
		return this->operator *(div);
	}
	else
	{
		// dividend includes zero
		if(right.inf==0 && right.R==0)
			return interval(interval::NEG_INF, 
							0,
							interval::POS_INF);
		if(this->sup<=0 && right.sup==0)
			return interval(this->sup/right.inf, 
							this->sup/right.inf,
							interval::POS_INF);
		if(this->sup<=0 && right.sup>0 && right.inf<0)
			return interval(interval::NEG_INF, 
							0,
							interval::POS_INF);
		if(this->sup<=0 && right.inf==0)
			return interval(interval::NEG_INF,
							this->sup/right.sup, 
							this->sup/right.sup);
		if(this->inf<0 && this->sup>0)
			return interval(interval::NEG_INF, 
							0,
							interval::POS_INF);
		if(this->inf>=0 && right.sup==0)
			return interval(interval::NEG_INF, 
							this->inf/right.inf,
							this->inf/right.inf);
		if(this->inf>=0 && right.sup>0 && right.inf<0)
			return interval(interval::NEG_INF, 
							0,
							interval::POS_INF);
		if(this->inf>=0 && right.inf==0)
			return interval(this->inf/right.sup, 
							this->inf/right.sup,
							interval::POS_INF);
		throw IVException("IVException: Divid by zero");
	}
*/
	if(a.INCLUDE(0))
		throw IVException("IVException: Divid by zero");
	else
	{
		a1b1 = inf/right.inf;
		a1b2 = inf/right.sup;
		a2b1 = sup/right.inf;
		a2b2 = sup/right.sup;

	if ( inf>=0 && sup>=0 && right.inf>0 && right.sup>0 )
		return interval(a1b2, a2b1);
	else if ( inf>=0 && sup>=0 && right.inf<0 && right.sup<0 )
		return interval(a2b2, a1b1);
	else if ( inf>=0 && sup<0 && right.inf>0 && right.sup>0 )
		return interval(a1b2, a2b2);
	else if ( inf>=0 && sup<0 && right.inf<0 && right.sup<0 )
		return interval(a2b1, a1b1);
	else if ( inf<0 && sup>=0 && right.inf>0 && right.sup>0 )
		return interval(a1b1, a2b1);
	else if ( inf<0 && sup>=0 && right.inf<0 && right.sup<0 )
		return interval(a2b2, a1b2);
	else if ( inf<0 && sup<0 && right.inf>0 && right.sup>0 )
		return interval(a1b1, a2b2);
	else if ( inf<0 && sup<0 && right.inf<0 && right.sup<0 )
		return interval(a2b1, a1b2);
	else
		throw IVException("IVException: Invalid division");
	}


}

interval interval::operator +(const double right)
{
	return interval(this->inf+right,this->sup+right);
}

interval interval::operator -(const double right)
{
	return interval(this->inf-right,this->sup-right);
}

interval interval::operator *(const double right)
{
	return interval(this->inf*right,this->sup*right);
}

interval interval::operator /(const double right)
{
	return interval(this->inf/right,this->sup/right);
}

interval& interval::operator=(const interval& right)
{
	if(this == &right) return *this;
	this->inf = right.inf;
	this->sup = right.sup;
	return *this;
}

interval& interval::operator+=(const interval& right)
{
	this->inf += right.inf;
	this->sup += right.sup;
	return *this;
}

interval& interval::operator-=(const interval& right)
{
	this->inf -= right.sup;
	this->sup -= right.inf;
	return *this;
}

interval& interval::operator*=(const interval& right)
{
	interval tmp = interval(this->inf,this->sup);
    tmp = tmp*right;
	this->inf = tmp.inf;
	this->sup = tmp.sup;    
	return *this;
}

interval& interval::operator/=(const interval& right)
{
	interval tmp = interval(this->inf,this->sup);
    tmp = tmp/right;
	this->inf = tmp.inf;
	this->sup = tmp.sup;    
	return *this;
}

interval& interval::operator=(const double right)
{
	this->inf = right;
	this->sup = right;
	return *this;
}

interval& interval::operator+=(const double right)
{
	this->inf += right;
	this->sup += right;
	return *this;
}

interval& interval::operator-=(const double right)
{
	this->inf -= right;
	this->sup -= right;
	return *this;
}

interval& interval::operator*=(const double right)
{
	this->inf *= right;
	this->sup *= right;    
	return *this;
}

interval& interval::operator/=(const double right)
{
	this->inf /= right;
	this->sup /= right;    
	return *this;
}


/*interval::operator int() const 
{
    return static_cast<int> this->Mid();
}*/

/*interval::operator double() const 
{
    return static_cast<double> this->Mid();
}
*/

interval interval::impr()
{
	if (inf >= sup)
		return interval(inf, sup);
	else
		return interval(sup, inf);
}

/*intervalItem interval::item()
{
		return intervalItem(inf, sup);
}*/

//////////////////////////////////////////////////////////////////////
//Friend mathmatical operations/functions
//////////////////////////////////////////////////////////////////////
interval operator+(double left, const interval& right)
{
	interval l(left);
	return l+right;
}

interval operator-(double left, const interval& right)
{
	interval l(left);
	return l-right;
}

interval operator*(double left, const interval& right)
{
	interval l(left);
	return l*right;
}

interval operator/(double left, const interval& right)
{
	interval l(left);
	return l/right;
}

interval operator-(const interval& right)
{
//	return interval(-right.sup, -right.nominal, -right.inf);
	return interval(-right.sup, -right.inf);
}

interval sin(const interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = sin(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return interval(-1,v_n,1);
		return interval(-1,1);
	else
	{
		double v_pi2 = PI/2;
		double v_3pi2 = 3*PI/2;
		double v_u = sin(p.sup);
		double v_l = sin(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(u>v_3pi2) u -= v_2pi;
		if(l>v_3pi2) l -= v_2pi;
		if(l<=u  &&  l<v_pi2  &&  u<v_pi2)
//			return interval(v_l, v_n, v_u);
			return interval(v_l, v_u);
		if(l<=u  &&  l<v_pi2  &&  u>=v_pi2)
			if(v_l <= v_u)
//				return interval(v_l, v_n, 1);
				return interval(v_l, 1);
			else 
//				return interval(v_u, v_n, 1);
				return interval(v_u, 1);
		if(l<=u  &&  l>=v_pi2)
//			return interval(v_u, v_n, v_l);
			return interval(v_u, v_l);
		if(l>u  &&  l<v_pi2  &&  u<v_pi2)
//			return interval(-1, v_n, 1);
			return interval(v_l, v_u);
		if(l>u  &&  l>=v_pi2  &&  u<v_pi2)
			if(v_l <= v_u)
//				return interval(-1, v_n, v_l);
				return interval(-1, v_l);
			else 
//				return interval(-1, v_n, v_u);
				return interval(-1, v_u);
		if(l>u  &&  u>=v_pi2)
//			return interval(-1, v_n, 1);
			return interval(v_u, v_l);
		throw new IVException("IVException: unknown at sin(interval)");
	}
}
/*
interval sin(const interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = sin(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return interval(-1,v_n,1);
		return interval(-1,1);
	else
	{
		double v_pi2 = PI/2;
		double v_3pi2 = 3*PI/2;
		double v_u = sin(p.sup);
		double v_l = sin(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(u>v_3pi2) u -= v_2pi;
		if(l>v_3pi2) l -= v_2pi;
		if(l<=u  &&  l<v_pi2  &&  u<v_pi2)
//			return interval(v_l, v_n, v_u);
			return interval(v_l, v_u);
		if(l<=u  &&  l<v_pi2  &&  u>=v_pi2)
			if(v_l <= v_u)
//				return interval(v_l, v_n, 1);
				return interval(v_l, 1);
			else 
//				return interval(v_u, v_n, 1);
				return interval(v_u, 1);
		if(l<=u  &&  l>=v_pi2)
//			return interval(v_u, v_n, v_l);
			return interval(v_u, v_l);
		if(l>u  &&  l<v_pi2  &&  u<v_pi2)
//			return interval(-1, v_n, 1);
			return interval(-1, 1);
		if(l>u  &&  l>=v_pi2  &&  u<v_pi2)
			if(v_l <= v_u)
//				return interval(-1, v_n, v_l);
				return interval(-1, v_l);
			else 
//				return interval(-1, v_n, v_u);
				return interval(-1, v_u);
		if(l>u  &&  u>=v_pi2)
//			return interval(-1, v_n, 1);
			return interval(-1, 1);
		throw new IVException("IVException: unknown at sin(interval)");
	}
}
*/

/*interval interval::power(double exp)
{
	double l, u;
	double a, b;

	if(exp>=0)
	{
		if(exp>1)
		{
			int i = exp/2;  
			if(i == exp/2.0) //if exp is even
			{
				if(inf>=0 && sup>=0)
				{
					l = pow(inf,exp);
					u = pow(sup,exp);
				}
				if(inf<0 && sup<0)
				{
					l = pow(sup,exp);
					u = pow(inf,exp);
				}
				if(inf<0 && sup>=0) 
				{
					a = pow(fabs(inf),exp);
					b = pow(fabs(sup),exp);
					l = 0.0;
					u = a>b?a:b;
				}
				if(inf>=0 && sup<0)
				{
					a = pow(fabs(inf),exp);
					b = pow(fabs(sup),exp);
					l = a>b?a:b;
					u = 0.0;
				}
			}
			else //if exp is odd
			{
				l = pow(inf,exp);
				u = pow(sup,exp);
			}
			return interval(l,u);
		}
		else
		{
			if(inf>=0 && sup>=0)
			{
				l = pow(inf,exp);
				u = pow(sup,exp);
				return interval(l,u);
			}
			else
				return interval(1,1);
		}


	}
	
	return 1/this->power(-exp);
}


interval interval::sqrt(interval v)
{
    return interval(sqrt(v.Inf()),sqrt(v.Sup()));
}*/

interval interval::dual()
{
//	return interval(sup, nominal, inf);
	return interval(sup, inf);
}

interval interval::prop()
{
	if (inf <= sup)
		return interval(inf, sup);
	else
		return interval(sup, inf);
}


interval cos(const interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = cos(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return interval(-1,v_n,1);
		return interval(-1,1);
	else
	{
		double v_u = cos(p.sup);
		double v_l = cos(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(l<=u  &&  l<PI  &&  u<PI)
//			return interval(v_u, v_n, v_l);
			return interval(v_u, v_l);
		if(l<=u  &&  l<PI  &&  u>=PI)
			if(v_l <= v_u)
//				return interval(-1, v_n, v_u);
				return interval(-1, v_u);
			else
//				return interval(-1, v_n, v_l);
				return interval(-1, v_l);
		if(l<=u  &&  l>=PI)
//			return interval(v_l, v_n, v_u);
			return interval(v_l, v_u);
		if(l>u  &&  l<PI  &&  u<PI)
//			return interval(-1, v_n, 1);
			return interval(v_u, v_l);
		if(l>u  &&  l>=PI  &&  u<PI)
			if(v_l <= v_u)
//				return interval(v_l, v_n, 1);
				return interval(v_l, 1);
			else 
//				return interval(v_u, v_n, 1);
				return interval(v_u, 1);
		if(l>u  &&  u>=PI)
//			return interval(-1, v_n, 1);
			return interval(v_l, v_u);
		throw new IVException("IVException: unknown at cos(interval)");
	}
}
/*
interval cos(const interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = cos(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return interval(-1,v_n,1);
		return interval(-1,1);
	else
	{
		double v_u = cos(p.sup);
		double v_l = cos(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(l<=u  &&  l<PI  &&  u<PI)
//			return interval(v_u, v_n, v_l);
			return interval(v_u, v_l);
		if(l<=u  &&  l<PI  &&  u>=PI)
			if(v_l <= v_u)
//				return interval(-1, v_n, v_u);
				return interval(-1, v_u);
			else
//				return interval(-1, v_n, v_l);
				return interval(-1, v_l);
		if(l<=u  &&  l>=PI)
//			return interval(v_l, v_n, v_u);
			return interval(v_l, v_u);
		if(l>u  &&  l<PI  &&  u<PI)
//			return interval(-1, v_n, 1);
			return interval(-1, 1);
		if(l>u  &&  l>=PI  &&  u<PI)
			if(v_l <= v_u)
//				return interval(v_l, v_n, 1);
				return interval(v_l, 1);
			else 
//				return interval(v_u, v_n, 1);
				return interval(v_u, 1);
		if(l>u  &&  u>=PI)
//			return interval(-1, v_n, 1);
			return interval(-1, 1);
		throw new IVException("IVException: unknown at cos(interval)");
	}
}
*/
//////////////////////////////////////////////////////////////////////
//Set operations
//////////////////////////////////////////////////////////////////////
interval interval::union_with(const interval& right)
{
	interval result = right;
/*	if(result.isEmpty())
	{
		result.setL(this->L);
		result.setNominal(this->nominal);
		result.setR(this->R);
	}
	else
*/	{
		if(this->inf < right.inf)
			result.setInf(this->inf);
		else
			result.setInf(right.inf);
		if(this->sup > right.sup)
			result.setSup(this->sup);
		else
			result.setSup(right.sup);
//		result.setNominal((this->nominal+right.nominal)/2);
	}
	return result;
}

interval interval::intersect_with(const interval& right)
{
	interval result = right;
/*	if(result.isEmpty())
	{
		result.setL(right.L);
		result.setNominal(right.nominal);
		result.setR(right.R);
	}
	else
*/	{
		if(this->inf > right.inf)
			result.setInf(this->inf);
		else
			result.setInf(right.inf);
		if(this->sup < right.sup)
			result.setSup(this->sup);
		else
			result.setSup(right.sup);
//		result.setNominal(result.Mid());
	}
	return result;
}


//////////////////////////////////////////////////////////////////////
// Relations
//////////////////////////////////////////////////////////////////////
int interval::operator ==(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if( inf==right.inf && 
//		nominal==right.nominal && 
		sup==right.sup)
		return 1;
	else
		return 0;
}

int interval::operator >=(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>=right.sup &&
		inf>=right.inf)
		return 1;
	else
		return 0;
}

int interval::operator >(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>right.sup &&
		inf>right.inf)
		return 1;
	else
		return 0;
}

int interval::operator <=(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<=right.sup &&
		inf<=right.inf)
		return 1;
	else
		return 0;
}

int interval::operator <(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<right.sup &&
		inf<right.inf)
		return 1;
	else
		return 0;
}

int interval::operator ==(const double right)
{
	if( inf==right && 
		sup==right)
		return 1;
	else
		return 0;
}

int interval::operator >=(const double right)
{
	if(sup>=right &&
		inf>=right)
		return 1;
	else
		return 0;
}

int interval::operator >(const double right)
{
	if(sup>right &&
		inf>right)
		return 1;
	else
		return 0;
}

int interval::operator <=(const double right)
{
	if(sup<=right &&
		inf<=right)
		return 1;
	else
		return 0;
}

int interval::operator <(const double right)
{
	if(sup<right &&
		inf<right)
		return 1;
	else
		return 0;
}

//////////////////////////////////////////////////////////////////////
// Comparison and relations
//////////////////////////////////////////////////////////////////////

int interval::EQ(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if( inf==right.inf && 
//		nominal==right.nominal && 
		sup==right.sup)
		return 1;
	else
		return 0;
}

int interval::GT_EQ(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>=right.sup &&
		inf>=right.inf)
		return 1;
	else
		return 0;
}

int interval::GT(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>right.sup &&
		inf>right.inf)
		return 1;
	else
		return 0;
}

int interval::LT_EQ(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<=right.sup &&
		inf<=right.inf)
		return 1;
	else
		return 0;
}

int interval::LT(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<right.sup &&
		inf<right.inf)
		return 1;
	else
		return 0;
}


int interval::S_GT_EQ(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(inf>=right.sup)
		return 1;
	else
		return 0;
}

int interval::S_GT(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(inf>right.sup)
		return 1;
	else
		return 0;
}

int interval::S_LT_EQ(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<=right.inf)
		return 1;
	else
		return 0;
}

int interval::S_LT(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<right.inf)
		return 1;
	else
		return 0;
}

int interval::INCLUDE(const interval &right)
{
/*	interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/
//	if( sup>=right.sup && inf<=right.inf  || sup<=right.sup && inf>=right.inf )
	if( sup>=right.sup && inf<=right.inf )
		return 1;
	else
		return 0;
}

