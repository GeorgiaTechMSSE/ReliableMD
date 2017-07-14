#include "imatrix.hpp"
#include "ivector.hpp"
#include "interval.hpp"  // predefined interval arithmetic
#include <iostream>
using namespace cxsc;
using namespace std;

int main()
{
  cout << "NOTE: Double-check Kaucher arithmetic in modified C-XSC library" << endl;
  interval x,y;
  x = interval(-3,7);
  y = interval(-2,-10);
  cout << " x = " << x << endl;
  cout << " y = " << y << endl;
  cout << "x*y = " << x*y << endl;
  cout << "diam(x) = " << diam(x) << endl;
  cout << "diam(y) = " << diam(y) << endl;
  cout << "radi(x) = " << radi(x) << endl;
  cout << "radi(y) = " << radi(y) << endl;
  interval z;
  z = 3;
  cout << " z = " << z << endl;

  // test interval array
  int n = 4;
  imatrix A(n,n);
  cout << "test interval matrix" << endl;
  for (int i = 1; i < n+1; i++)
    for (int j = 1; j < n+1; j++) {
      A[i][j] = interval(j,i);
      cout << "A[" << i << "][" << j << "] = " << A[i][j] << endl;
    }

  cout << endl;
  cout << endl;

  cout << "test interval vector" << endl;

  ivector Z(n);
  for (int i = 1; i < n+1; i++) {
    Z[i] = i;
    cout << "Z[" << i << "] = " << Z[i] << endl;
  }

  cout << endl;
  cout << endl;

  ivector AtimesZ(n);
  AtimesZ = A*Z;
  for (int i = 1; i < n+1; i++ )
    cout << "AtimesZ[" << i << "] =" << AtimesZ[i] << endl;

  return 0;
}