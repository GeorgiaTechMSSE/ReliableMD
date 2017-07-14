// #include "example.hpp"
#include "real.hpp"
#include "imatrix.hpp"
#include "ivector.hpp"
#include "interval.hpp"  // predefined interval arithmetic
#include "rvector.hpp"
#include "l_interval.hpp"
#include <iostream>
#include <ctime>
using namespace cxsc;
using namespace std;

int main()
{
	interval x = interval(1,2);
	cout << x*x*x << endl;
	cout << interval(1.0)/interval(1,2) << endl;
}