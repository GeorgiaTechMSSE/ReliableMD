//============================================================================
//
// 	      ParLinSys - A Solver for Parametric Linear Systems
//
//	    Supplement to the C-XSC Toolbox for Verified Computing
//
// Author: Michael Zimmer, based on a previous version by Evgenija Popova
//
// This program is free software for non-commercial use.
//
// For details on theory see the papers:
//
// S. Rump: Verification methods for dense and sparse systems of equations,
// in: Topics in Validated Computations, ed. J. Herzberger, Elsevier Science
// B.V., 1994,	pp. 63-135.
//
// E. Popova, W. Kraemer: Parametric Fixed-Point Iteration Implemented in
// C-XSC. Universitaet Wuppertal, Preprint BUW - WRSWT, 2003.
//
// M. Zimmer, W. Kraemer, E.Popova: Solvers for the verified solution of 
// parametric linear systems, Computing: Volume 94, Issue 2 (2012), Page 109-123,
// Springer.
//
// For versions using  the alternative algorithm:
// A. Neumaier, A. Pownuk: Linear systems wit hlarge uncertainties, with 
// applications to truss structures. Reliable Computing, 13(149-172), 2007.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================


/*
 *  This header contains the helper classe CoeffMatrix, CoeffIMatrix, CoeffCMatrix
 *  and CoeffCIMatrix. These are wrapper classes for either sparse or full matrices
 *  of the according data type and are used for the coefficient matrices in the 
 *  parametric solver. See example files for usage of these classes.
 */

#include <iostream>
#include <scimatrix.hpp>

namespace cxsc {
  
static const rmatrix& Re(const rmatrix& A) {
  return A;
}

static rmatrix Im(const rmatrix& A) {
  rmatrix R(ColLen(A),RowLen(A));
  R = 0.0;
  return R;
}

static const srmatrix& Re(const srmatrix& A) {
  return A;
}

static srmatrix Im(const srmatrix& A) {
  rmatrix R(ColLen(A),RowLen(A));
  return R;
}


class CoeffMatrix;
class CoeffCMatrix;
class CoeffCIMatrix;
  
class CoeffIMatrix {
  private:
    cxsc::imatrix F;
    cxsc::simatrix S;
    bool fullmatrix;
    
    
  public:
    CoeffIMatrix() : fullmatrix(true) {}
    CoeffIMatrix(const cxsc::imatrix& A) : fullmatrix(true), F(A) {}
    CoeffIMatrix(const cxsc::simatrix& A) : fullmatrix(false), S(A) {}
    CoeffIMatrix(const cxsc::CoeffMatrix& A);
    
    bool isFull() const { return fullmatrix; }
    cxsc::imatrix& full() { return F; }
    cxsc::simatrix& sparse() { return S; }
    const cxsc::imatrix& full() const { return F; }
    const cxsc::simatrix& sparse() const { return S; }

    friend cxsc::imatrix& operator+= (cxsc::imatrix& m1, const CoeffIMatrix& m2);
    friend cxsc::imatrix& operator-= (cxsc::imatrix& m1, const CoeffIMatrix& m2);
    friend cxsc::imatrix operator* (const cxsc::rmatrix& m1, const CoeffIMatrix& m2);
    
    friend class CoeffCIMatrix;
};


inline cxsc::imatrix& operator+= (cxsc::imatrix& m1, const CoeffIMatrix& m2) {
  if(m2.fullmatrix)
    return m1+=m2.F;
  else
    return m1+=m2.S;
}

inline cxsc::imatrix& operator-= (cxsc::imatrix& m1, const CoeffIMatrix& m2) {
  if(m2.fullmatrix)
    return m1-=m2.F;
  else
    return m1-=m2.S;
}

inline cxsc::imatrix operator* (const cxsc::rmatrix& m1, const CoeffIMatrix& m2) {
  if(m2.fullmatrix)
    return m1 * m2.F;
  else
    return m1 * m2.S;
}
    
    

class CoeffMatrix {
  private:
    cxsc::rmatrix F;
    cxsc::srmatrix S;
    bool fullmatrix;
    
    
  public:
    CoeffMatrix() : fullmatrix(true) {}
    CoeffMatrix(const cxsc::rmatrix& A) : fullmatrix(true), F(A) {}
    CoeffMatrix(const cxsc::srmatrix& A) : fullmatrix(false), S(A) {}

    bool isFull() const { return fullmatrix; }
    cxsc::rmatrix& full() { return F; }
    cxsc::srmatrix& sparse() { return S; }
    const cxsc::rmatrix& full() const { return F; }
    const cxsc::srmatrix& sparse() const { return S; }

    friend std::ostream& operator<<(std::ostream& os, const CoeffMatrix& c);
    friend CoeffMatrix operator+ (const CoeffMatrix& m1, const CoeffMatrix& m2);
    friend CoeffMatrix operator- (const CoeffMatrix& m1, const CoeffMatrix& m2);
    friend CoeffMatrix operator* (const CoeffMatrix& m1, const CoeffMatrix& m2);
    friend CoeffMatrix operator* (const cxsc::rmatrix& m1, const CoeffMatrix& m2);
    friend cxsc::rvector operator* (const CoeffMatrix& m, const cxsc::rvector& v);
    friend cxsc::cvector operator* (const CoeffMatrix& m, const cxsc::cvector& v);
    friend cxsc::ivector operator* (const CoeffMatrix& m, const cxsc::ivector& v);
    friend cxsc::civector operator* (const CoeffMatrix& m, const cxsc::civector& v);
    friend CoeffIMatrix operator* (const cxsc::interval& i, const CoeffMatrix& m);
    friend CoeffIMatrix operator* (const CoeffMatrix& m, const cxsc::interval& i);
    friend CoeffCIMatrix operator* (const CoeffMatrix& m, const cxsc::cinterval& i);
    friend CoeffMatrix operator* (const cxsc::real& i, const CoeffMatrix& m);
    friend CoeffCMatrix operator* (const cxsc::complex& i, const CoeffMatrix& m);
    friend CoeffCIMatrix operator* (const cxsc::cinterval& i, const CoeffMatrix& m);
    friend cxsc::rmatrix& operator+= (cxsc::rmatrix& m1, const CoeffMatrix& m2);    
//     friend cxsc::rmatrix& operator=(cxsc::rmatrix& i, const CoeffMatrix& m);
//     friend cxsc::imatrix& operator=(cxsc::imatrix& i, const CoeffMatrix& m);    
    friend cxsc::srmatrix Id(CoeffMatrix& m);
    friend int RowLen(CoeffMatrix& m);
    friend int ColLen(CoeffMatrix& m);   
    
    friend class CoeffIMatrix;
    friend class CoeffCIMatrix;
};

    inline std::ostream& operator<<(std::ostream& os, const CoeffMatrix& c) {
        if(c.isFull()) {
            os << c.full();
        } else {
            os << c.sparse();
        }
        return os;
    }


inline CoeffMatrix operator+ (const CoeffMatrix& m1, const CoeffMatrix& m2) {
  if(m1.fullmatrix && m2.fullmatrix) 
    return CoeffMatrix(m1.F + m2.F);
  else if(m1.fullmatrix && !m2.fullmatrix)
    return CoeffMatrix(m1.F + m2.S);
  else if(!m1.fullmatrix && m2.fullmatrix)
    return CoeffMatrix(m1.S + m2.F);
  else
    return CoeffMatrix(m1.S + m2.S);    
}

inline CoeffMatrix operator- (const CoeffMatrix& m1, const CoeffMatrix& m2) {
  if(m1.fullmatrix && m2.fullmatrix) 
    return CoeffMatrix(m1.F - m2.F);
  else if(m1.fullmatrix && !m2.fullmatrix)
    return CoeffMatrix(m1.F - m2.S);
  else if(!m1.fullmatrix && m2.fullmatrix)
    return CoeffMatrix(m1.S - m2.F);
  else
    return CoeffMatrix(m1.S - m2.S);    
}

inline CoeffMatrix operator* (const CoeffMatrix& m1, const CoeffMatrix& m2) {
  if(m1.fullmatrix && m2.fullmatrix) 
    return CoeffMatrix(m1.F * m2.F);
  else if(m1.fullmatrix && !m2.fullmatrix)
    return CoeffMatrix(m1.F * m2.S);
  else if(!m1.fullmatrix && m2.fullmatrix)
    return CoeffMatrix(m1.S * m2.F);
  else
    return CoeffMatrix(m1.S * m2.S);    
}

inline CoeffMatrix operator* (const cxsc::rmatrix& m1, const CoeffMatrix& m2){
  if(m2.fullmatrix)
    return CoeffMatrix(m1*m2.F);
  else
    return CoeffMatrix(m1*m2.S);
}

inline cxsc::rvector operator* (const CoeffMatrix& m, const cxsc::rvector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline cxsc::cvector operator* (const CoeffMatrix& m, const cxsc::cvector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline cxsc::ivector operator* (const CoeffMatrix& m, const cxsc::ivector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline CoeffIMatrix operator* (const cxsc::interval& i, const CoeffMatrix& m) {
  if(m.fullmatrix)
    return CoeffIMatrix(i*m.F);
  else
    return CoeffIMatrix(i*m.S);
}

inline CoeffIMatrix operator* (const CoeffMatrix& m, const cxsc::interval& i) {
  if(m.fullmatrix)
    return CoeffIMatrix(i*m.F);
  else
    return CoeffIMatrix(i*m.S);
}

inline CoeffMatrix operator* (const cxsc::real& i, const CoeffMatrix& m) {
  if(m.fullmatrix)
    return CoeffMatrix(i*m.F);
  else
    return CoeffMatrix(i*m.S);
}

inline cxsc::rmatrix& operator+= (cxsc::rmatrix& m1, const CoeffMatrix& m2) {
  if(m2.fullmatrix)
    return m1 += m2.F;
  else
    return m1 += m2.S;
}

// cxsc::rmatrix& operator=(cxsc::rmatrix& i, const CoeffMatrix& m) {
//   if(m.fullmatrix)
//     return i = m.F;
//   else
//     return i = m.S;
// }
// 
// cxsc::imatrix& operator=(cxsc::imatrix& i, const CoeffMatrix& m) {
//   if(m.fullmatrix)
//     return i = m.F;
//   else
//     return i = m.S;
// }

inline cxsc::srmatrix Id(CoeffMatrix& m) {
  if(m.fullmatrix)
    return Id(m.F);
  else
    return Id(m.S);
}

inline int RowLen(CoeffMatrix& m) {
  if(m.fullmatrix)
    return RowLen(m.F);
  else
    return RowLen(m.S);
  
}

inline int ColLen(CoeffMatrix& m) {
  if(m.fullmatrix)
    return ColLen(m.F);
  else
    return ColLen(m.S);
}

inline CoeffIMatrix::CoeffIMatrix(const cxsc::CoeffMatrix& A) : fullmatrix(A.fullmatrix), F(A.F), S(A.S) { }


class CoeffCMatrix;
  
class CoeffCIMatrix {
  private:
    cxsc::cimatrix F;
    cxsc::scimatrix S;
    bool fullmatrix;
    
    
  public:
    CoeffCIMatrix() : fullmatrix(true) {}
    CoeffCIMatrix(const cxsc::cimatrix& A) : fullmatrix(true), F(A) {}
    CoeffCIMatrix(const cxsc::scimatrix& A) : fullmatrix(false), S(A) {}
    CoeffCIMatrix(const cxsc::imatrix& A) : fullmatrix(true), F(A) {}
    CoeffCIMatrix(const cxsc::simatrix& A) : fullmatrix(false), S(A) {}
    CoeffCIMatrix(const cxsc::CoeffCMatrix& A);
    CoeffCIMatrix(const cxsc::CoeffMatrix& A);
    CoeffCIMatrix(const cxsc::CoeffIMatrix& A);
    
    bool isFull() const { return fullmatrix; }
    cxsc::cimatrix& full() { return F; }
    cxsc::scimatrix& sparse() { return S; }
    const cxsc::cimatrix& full() const { return F; }
    const cxsc::scimatrix& sparse() const { return S; }


    friend cxsc::cimatrix& operator+= (cxsc::imatrix& m1, const CoeffCIMatrix& m2);
    friend cxsc::cimatrix operator* (const cxsc::rmatrix& m1, const CoeffCIMatrix& m2);
    friend cxsc::cimatrix& operator+= (cxsc::cimatrix& m1, const CoeffCIMatrix& m2);
    friend cxsc::cimatrix& operator-= (cxsc::cimatrix& m1, const CoeffCIMatrix& m2);
    friend cxsc::cimatrix operator* (const cxsc::cmatrix& m1, const CoeffCIMatrix& m2);
};

inline cxsc::cimatrix& operator+= (cxsc::imatrix& m1, const CoeffCIMatrix& m2) {
  if(m2.fullmatrix)
    return m1+=m2.F;
  else
    return m1+=m2.S;
}

inline cxsc::cimatrix& operator-= (cxsc::cimatrix& m1, const CoeffCIMatrix& m2) {
  if(m2.fullmatrix)
    return m1-=m2.F;
  else
    return m1-=m2.S;
}

inline cxsc::cimatrix operator* (const cxsc::rmatrix& m1, const CoeffCIMatrix& m2) {
  if(m2.fullmatrix)
    return m1 * m2.F;
  else
    return m1 * m2.S;
}

inline cxsc::cimatrix& operator+= (cxsc::cimatrix& m1, const CoeffCIMatrix& m2) {
  if(m2.fullmatrix)
    return m1+=m2.F;
  else
    return m1+=m2.S;
}

inline cxsc::cimatrix operator* (const cxsc::cmatrix& m1, const CoeffCIMatrix& m2) {
  if(m2.fullmatrix)
    return m1 * m2.F;
  else
    return m1 * m2.S;
}

    

class CoeffCMatrix {
  private:
    cxsc::cmatrix F;
    cxsc::scmatrix S;
    bool fullmatrix;
    
    
  public:
    CoeffCMatrix() : fullmatrix(true) {}
    CoeffCMatrix(const cxsc::rmatrix& A) : fullmatrix(true), F(A) {}
    CoeffCMatrix(const cxsc::srmatrix& A) : fullmatrix(false), S(A) {}
    CoeffCMatrix(const cxsc::cmatrix& A) : fullmatrix(true), F(A) {}
    CoeffCMatrix(const cxsc::scmatrix& A) : fullmatrix(false), S(A) {}

    bool isFull() const { return fullmatrix; }
    cxsc::cmatrix& full() { return F; }
    cxsc::scmatrix& sparse() { return S; }
    const cxsc::cmatrix& full() const { return F; }
    const cxsc::scmatrix& sparse() const { return S; }
    

    friend std::ostream& operator<<(std::ostream& os, const CoeffCMatrix& c);
    friend CoeffCMatrix operator+ (const CoeffCMatrix& m1, const CoeffCMatrix& m2);
    friend CoeffCMatrix operator- (const CoeffCMatrix& m1, const CoeffCMatrix& m2);
    friend CoeffCMatrix operator* (const CoeffCMatrix& m1, const CoeffCMatrix& m2);
    friend CoeffCMatrix operator* (const cxsc::rmatrix& m1, const CoeffCMatrix& m2);
    friend CoeffCMatrix operator* (const cxsc::cmatrix& m1, const CoeffCMatrix& m2);
    friend cxsc::cvector operator* (const CoeffCMatrix& m, const cxsc::rvector& v);
    friend cxsc::civector operator* (const CoeffCMatrix& m, const cxsc::ivector& v);
    friend cxsc::cvector operator* (const CoeffCMatrix& m, const cxsc::cvector& v);
    friend cxsc::civector operator* (const CoeffCMatrix& m, const cxsc::civector& v);
    friend CoeffCMatrix operator* (const cxsc::real& i, const CoeffCMatrix& m);
    friend CoeffCIMatrix operator* (const cxsc::cinterval& i, const CoeffCMatrix& m);
        friend CoeffCIMatrix operator* (const CoeffCMatrix& m, const cxsc::cinterval& i);
    friend CoeffCMatrix operator* (const cxsc::complex& i, const CoeffCMatrix& m);
    friend cxsc::cmatrix& operator+= (cxsc::cmatrix& m1, const CoeffCMatrix& m2);    
//     friend cxsc::rmatrix& operator=(cxsc::rmatrix& i, const CoeffCMatrix& m);
//     friend cxsc::imatrix& operator=(cxsc::imatrix& i, const CoeffCMatrix& m);    
    friend cxsc::srmatrix Id(CoeffCMatrix& m);
    friend int RowLen(CoeffCMatrix& m);
    friend int ColLen(CoeffCMatrix& m);   
    
    friend class CoeffCIMatrix;
};

    inline std::ostream& operator<<(std::ostream& os, const CoeffCMatrix& c) {
        if(c.isFull()) {
            os << c.full();
        } else {
            os << c.sparse();
        }
        return os;
    }


inline CoeffCMatrix operator+ (const CoeffCMatrix& m1, const CoeffCMatrix& m2) {
  if(m1.fullmatrix && m2.fullmatrix) 
    return CoeffCMatrix(m1.F + m2.F);
  else if(m1.fullmatrix && !m2.fullmatrix)
    return CoeffCMatrix(m1.F + m2.S);
  else if(!m1.fullmatrix && m2.fullmatrix)
    return CoeffCMatrix(m1.S + m2.F);
  else
    return CoeffCMatrix(m1.S + m2.S);    
}

inline CoeffCMatrix operator- (const CoeffCMatrix& m1, const CoeffCMatrix& m2) {
  if(m1.fullmatrix && m2.fullmatrix) 
    return CoeffCMatrix(m1.F - m2.F);
  else if(m1.fullmatrix && !m2.fullmatrix)
    return CoeffCMatrix(m1.F - m2.S);
  else if(!m1.fullmatrix && m2.fullmatrix)
    return CoeffCMatrix(m1.S - m2.F);
  else
    return CoeffCMatrix(m1.S - m2.S);    
}

inline CoeffCMatrix operator* (const CoeffCMatrix& m1, const CoeffCMatrix& m2) {
  if(m1.fullmatrix && m2.fullmatrix) 
    return CoeffCMatrix(m1.F * m2.F);
  else if(m1.fullmatrix && !m2.fullmatrix)
    return CoeffCMatrix(m1.F * m2.S);
  else if(!m1.fullmatrix && m2.fullmatrix)
    return CoeffCMatrix(m1.S * m2.F);
  else
    return CoeffCMatrix(m1.S * m2.S);    
}

inline CoeffCMatrix operator* (const cxsc::rmatrix& m1, const CoeffCMatrix& m2){
  if(m2.fullmatrix)
    return CoeffCMatrix(m1*m2.F);
  else
    return CoeffCMatrix(m1*m2.S);
}

inline CoeffCMatrix operator* (const cxsc::cmatrix& m1, const CoeffCMatrix& m2){
  if(m2.fullmatrix)
    return CoeffCMatrix(m1*m2.F);
  else
    return CoeffCMatrix(m1*m2.S);
}

inline cxsc::cvector operator* (const CoeffCMatrix& m, const cxsc::rvector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline cxsc::civector operator* (const CoeffCMatrix& m, const cxsc::ivector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline cxsc::cvector operator* (const CoeffCMatrix& m, const cxsc::cvector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline cxsc::civector operator* (const CoeffCMatrix& m, const cxsc::civector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline CoeffCMatrix operator* (const cxsc::real& i, const CoeffCMatrix& m) {
  if(m.fullmatrix)
    return CoeffCMatrix(i*m.F);
  else
    return CoeffCMatrix(i*m.S);
}

inline CoeffCIMatrix operator* (const cxsc::cinterval& i, const CoeffCMatrix& m) {
  if(m.fullmatrix)
    return CoeffCIMatrix(i*m.F);
  else
    return CoeffCIMatrix(i*m.S);
}

inline CoeffCIMatrix operator* (const CoeffCMatrix& m, const cxsc::cinterval& i) {
  if(m.fullmatrix)
    return CoeffCIMatrix(i*m.F);
  else
    return CoeffCIMatrix(i*m.S);
}

inline CoeffCMatrix operator* (const cxsc::complex& i, const CoeffCMatrix& m) {
  if(m.fullmatrix)
    return CoeffCMatrix(i*m.F);
  else
    return CoeffCMatrix(i*m.S);
}

inline cxsc::cmatrix& operator+= (cxsc::cmatrix& m1, const CoeffCMatrix& m2) {
  if(m2.fullmatrix)
    return m1 += m2.F;
  else
    return m1 += m2.S;
}

// cxsc::rmatrix& operator=(cxsc::rmatrix& i, const CoeffCMatrix& m) {
//   if(m.fullmatrix)
//     return i = m.F;
//   else
//     return i = m.S;
// }
// 
// cxsc::imatrix& operator=(cxsc::imatrix& i, const CoeffCMatrix& m) {
//   if(m.fullmatrix)
//     return i = m.F;
//   else
//     return i = m.S;
// }

inline cxsc::srmatrix Id(CoeffCMatrix& m) {
  if(m.fullmatrix)
    return Id(Re(m.F));
  else
    return Id(Re(m.S));
}

inline int RowLen(CoeffCMatrix& m) {
  if(m.fullmatrix)
    return RowLen(m.F);
  else
    return RowLen(m.S);
  
}

inline int ColLen(CoeffCMatrix& m) {
  if(m.fullmatrix)
    return ColLen(m.F);
  else
    return ColLen(m.S);
}

inline cxsc::civector operator* (const CoeffMatrix& m, const cxsc::civector& v) {
  if(m.fullmatrix)
    return m.F*v;
  else
    return m.S*v;
}

inline CoeffCIMatrix operator* (const cxsc::cinterval& i, const CoeffMatrix& m) {
  if(m.fullmatrix)
    return CoeffCIMatrix(i*m.F);
  else
    return CoeffCIMatrix(i*m.S);
}

inline CoeffCMatrix operator* (const cxsc::complex& i, const CoeffMatrix& m) {
  if(m.fullmatrix)
    return CoeffCMatrix(i*m.F);
  else
    return CoeffCMatrix(i*m.S);
}

inline CoeffCIMatrix operator* (const CoeffMatrix& m, const cxsc::cinterval& i) {
  if(m.fullmatrix)
    return CoeffCIMatrix(i*m.F);
  else
    return CoeffCIMatrix(i*m.S);
}


inline CoeffCIMatrix::CoeffCIMatrix(const cxsc::CoeffCMatrix& A) : fullmatrix(A.fullmatrix), F(A.F), S(A.S) { }
inline CoeffCIMatrix::CoeffCIMatrix(const cxsc::CoeffMatrix& A) : fullmatrix(A.fullmatrix), F(A.F), S(A.S) { }
inline CoeffCIMatrix::CoeffCIMatrix(const cxsc::CoeffIMatrix& A) : fullmatrix(A.fullmatrix), F(A.F), S(A.S) { }

  
} //namespace cxsc