#ifndef HAVE_MATRIX_HPP
#define HAVE_MATRIX_HPP

#include<algorithm>
#include<array>
#include<cmath>
#include<complex>
#include<memory>

#ifdef VERBOSE
#  define DOIFVERBOSE( S ) S
#else
#  define DOIFVERBOSE( S ) // S 
#endif

#ifdef MOVE
#  define DOIFMOVE( S ) S
#else
#  define DOIFMOVE( S ) // S 
#endif

#ifdef PRINTRES
#  define DOIFPRINTRES( S ) S
#else
#  define DOIFPRINTRES( S ) // S 
#endif

// Operator+ overloading to add complex and real
template<class Scalar1, class Scalar2>
auto operator+ (const Scalar1& m1, const std::complex<Scalar2>& m2) {
  using return_type = std::complex<decltype (Scalar1{} + Scalar2{})>;
  return_type res(m1+m2.real (), m2.imag ());
  return res;
}

template<class Scalar1, class Scalar2>
auto operator+ (const std::complex<Scalar2>& m2, const Scalar1& m1) {
  using return_type = std::complex<decltype (Scalar1{} + Scalar2{})>;
  return_type res(m1+m2.real (), m2.imag ());
  return res;
}

// Operator* overloading to multiply complex and real
template<class Scalar1, class Scalar2>
auto operator* (const Scalar1& m1, const std::complex<Scalar2>& m2) {
  using return_type = std::complex<decltype (Scalar1{} * Scalar2{})>;
  return_type res(m1 * m2.real (), m1 * m2.imag ());
  return res;
}

template<class Scalar1, class Scalar2>
auto operator* (const std::complex<Scalar2>& m2, const Scalar1& m1) {
  using return_type = std::complex<decltype (Scalar1{} + Scalar2{})>;
  return_type res(m1 * m2.real (), m1 * m2.imag ());
  return res;
}

// Operator/ overloading to divide complex by real
template<class Scalar1, class Scalar2>
auto operator/ (const std::complex<Scalar1>& m1, const Scalar2& m2) {
  using return_type = std::complex<decltype (Scalar1{} / Scalar2{})>;
  return_type res(m1.real () / m2, m1.imag () / m2);
  return res;
}

// Operator/ overloading to divide real by complex
template<class Scalar1, class Scalar2>
auto operator/ (const Scalar1& m1, const std::complex<Scalar2>& m2) {
  using return_type = std::complex<decltype (Scalar1{} / Scalar2{})>;
  return_type res = return_type(m1, 0) / return_type(m2.real (), m2.imag ());
  return res;
}

template<class Scalar, std::size_t Rows, std::size_t Cols>
class Matrix {

private :
  
  //std::unique_ptr<Scalar[]> data;
  std::vector<Scalar> data;
  
public :

  template<class S, std::size_t R>
  using ColumnVector = Matrix<S, R, 1>;

  Matrix() : data(Rows*Cols, static_cast<Scalar> (0)) {
    DOIFVERBOSE ( std::cout << "calling default ctor" << " data.data () == " << data.data () << std::endl; )
      //std::fill_n (data.get (), Rows*Cols, Scalar{0});
  };
  
  ~Matrix() {
    DOIFVERBOSE ( std::cout << "calling default dctor" << " data.data () == " << data.data () << std::endl; )
  };
  
  Scalar* get_data () {
    return data.data ();
  };

  const Scalar* get_data () const {
    return data.data ();
  };


  DOIFMOVE (
  Matrix (Matrix&& inm) {
    DOIFVERBOSE ( std::cout << "calling move ctor" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; );
     this->data.swap (inm.data);
    (std::vector<Scalar>{}).swap (inm.data);
    DOIFVERBOSE ( std::cout << "move ctor done" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; );
  }
            );
  
  
  // Matrix (const Matrix&) = default; // does not work as long as data is a unique_ptr!
  Matrix (const Matrix& inm) : data(Rows*Cols, static_cast<Scalar> (0)) {
    DOIFVERBOSE ( std::cout << "calling copy ctor" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; )
    for (std::size_t i = 0; i < Rows; ++i) 
      for (std::size_t j = 0; j < Cols; ++j) 
        (*this)(i, j) = inm (i, j);
  }
  
  // copy constructor with data conversion
  template<class InScalar>
  Matrix (const Matrix<InScalar, Rows, Cols>& inm) : data(Rows*Cols, static_cast<Scalar> (0)) {
    DOIFVERBOSE ( std::cout << "calling copy ctor" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; )
    for (std::size_t i = 0; i < Rows; ++i)
      for (std::size_t j = 0; j < Cols; ++j) 
        (*this)(i, j) = static_cast<Scalar> (inm (i, j));
  }

  // copying complex matrix to real matrix silently discards imaginary part
  template<class InScalar>
  Matrix (const Matrix<std::complex<InScalar>, Rows, Cols>& inm) : data(Rows*Cols, static_cast<Scalar> (0)) {
    DOIFVERBOSE ( std::cout << "calling copy ctor" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; )
    for (std::size_t i = 0; i < Rows; ++i) 
      for (std::size_t j = 0; j < Cols; ++j) 
        (*this)(i, j) = (inm (i, j)).real ();
  }
  
  // Matrix& operator= (const Matrix&) = default; // does not work as long as data is a unique_ptr!
  Matrix& operator= (const Matrix& inm) {
    DOIFVERBOSE ( std::cout << "calling operator= (const Matrix& inm)" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; )
    for (std::size_t i = 0; i < Rows; ++i) 
      for (std::size_t j = 0; j < Cols; ++j) 
        (*this)(i, j) = inm (i, j);
    return (*this);
  }


  template<class InScalar>
  Matrix& operator= (const Matrix<InScalar, Rows, Cols>& inm) {
    DOIFVERBOSE ( std::cout << "calling operator= (const Matrix<InScalar, Rows, Cols>& inm)" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; )
    for (std::size_t i = 0; i < Rows; ++i) 
      for (std::size_t j = 0; j < Cols; ++j) 
        (*this)(i, j) = static_cast<Scalar> (inm (i, j));
    return (*this);
  }

  // assigning complex matrix to real matrix silently discards imaginary part
  template<class InScalar>
  Matrix& operator= (const Matrix<std::complex<InScalar>, Rows, Cols>& inm) {
    DOIFVERBOSE ( std::cout << "calling operator= (const Matrix<std::complex<InScalar>, Rows, Cols>& inm)" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; )
    for (std::size_t i = 0; i < Rows; ++i) {
      for (std::size_t j = 0; j < Cols; ++j) {
        (*this)(i, j) = (inm (i, j)).real ();
      }
    }
    return (*this);
  }

  // If not doing data conversion we can also move on assignment
  DOIFMOVE (
  Matrix& operator= (Matrix&& inm) {
    DOIFVERBOSE ( std::cout << "calling operator= (Matrix&& inm)" << " data.data () == " << data.data () << " inm.data.data () == " << inm.data.data () << std::endl; ) ;
    this->data.swap (inm.data);
    (std::vector<Scalar>{}).swap (inm.data);
    return (*this);
  }
            );
  
  Scalar&
  operator() (std::size_t i, std::size_t j = 0) {
    return data[i+j*Rows];
  }

  const Scalar&
  operator() (std::size_t i, std::size_t j = 0) const {
    return data[i+j*Rows];
  }

  Matrix<Scalar, Cols, Rows>
  transpose () const {
    Matrix<Scalar, Cols, Rows> mtmp{};
    for (std::size_t i = 0; i < Rows; ++i)
      for (std::size_t j = 0; j < Cols; ++j)
        mtmp(j, i) = (*this) (i, j);
    return mtmp;
  }

  Matrix
  operator- () const {
    DOIFVERBOSE ( std::cout << "calling unary operator-" << std::endl; )
    Matrix mtmp{};
    for (std::size_t i = 0; i < Rows; ++i)
      for (std::size_t j = 0; j < Cols; ++j)
        mtmp(i, j) =  -(*this) (i, j);
    return mtmp;
  }
  
  void
  print () const {
    for (std::size_t i = 0; i < Rows; ++i) {
      for (std::size_t j = 0; j < Cols; ++j) 
        std::cout << (*this)(i,j) << " ";
      std::cout << std::endl;
    }
  }

  //  ColumnVector<Scalar, Rows>
  // solve (const ColumnVector<Scalar, Rows>& b);

  template<class Scalar2> 
  auto
  solve (const ColumnVector<Scalar2, Rows>& b) {
    DOIFVERBOSE ( std::cout << "calling method solve()" << std::endl; )
      auto x = solve_helper (*this, b);
    std::cout << "method solve done" << std::endl;
    return x;
  }
  
};



template<class Scalar1, class Scalar2, std::size_t Rows, std::size_t Cols>
auto operator+ (const Matrix<Scalar1, Rows, Cols>& m1, const Matrix<Scalar2, Rows, Cols>& m2) {
  DOIFVERBOSE ( std::cout << "calling operator+" << std::endl; )
  using return_type = Matrix<decltype (m1(0,0) + m2(0,0)), Rows, Cols>;
  return_type mtmp{};
  for (std::size_t i = 0; i < Rows; ++i)
    for (std::size_t j = 0; j < Cols; ++j) {
      mtmp(i,j) = m1(i,j) + m2(i, j);
    }

  std::cout << "operator+ done" << std::endl;
  return mtmp;
}

// Specialization for matrices of the same type
template<class Scalar, std::size_t Rows, std::size_t Cols>
Matrix<Scalar, Rows, Cols> operator+ (const Matrix<Scalar, Rows, Cols>& m1, const Matrix<Scalar, Rows, Cols>& m2) {
  DOIFVERBOSE ( std::cout << "calling operator+" << std::endl; ) ;
  Matrix<Scalar, Rows, Cols> mtmp{};
  for (std::size_t i = 0; i < Rows; ++i)
    for (std::size_t j = 0; j < Cols; ++j) {
      mtmp(i,j) = m1(i,j) + m2(i, j);
    }

  std::cout << "operator+ done" << std::endl;
  return mtmp;
}

template<class Scalar1, class Scalar2, std::size_t Rows, std::size_t Cols, std::size_t N>
auto operator* (const Matrix<Scalar1, Rows, N>& m1, const Matrix<Scalar2, N, Cols>& m2) {
  DOIFVERBOSE ( std::cout << "calling operator*" << std::endl; )
  using return_type = Matrix<decltype (m1(0,0) * m2(0,0)), Rows, Cols>;

  return_type mtmp{};
   
  for (std::size_t     i = 0; i < Rows; ++i)
    for (std::size_t   j = 0; j < Cols; ++j)
      for (std::size_t k = 0; k < N;    ++k) 
        mtmp(i,j) += m1(i, k) * m2(k, j);

  DOIFVERBOSE ( std::cout << "operator* done" << std::endl; )
  return mtmp;
}

template<class Scalar1, class Scalar2, std::size_t Rows, std::size_t Cols>
auto operator- (const Matrix<Scalar1, Rows, Cols>& m1, const Matrix<Scalar2, Rows, Cols>& m2) {     
  DOIFVERBOSE ( std::cout << "calling operator-" << std::endl; )
  auto m = m1 + (- m2);
  DOIFVERBOSE ( std::cout << "operator- done" << std::endl; )
  return m;
}

// Specialization for matrices of the same type
template<class Scalar, std::size_t Rows, std::size_t Cols>
Matrix<Scalar, Rows, Cols> operator- (const Matrix<Scalar, Rows, Cols>& m1, const Matrix<Scalar, Rows, Cols>& m2) {     
  DOIFVERBOSE ( std::cout << "calling operator-" << std::endl; )
  auto m = m1 + (- m2);
  DOIFVERBOSE ( std::cout << "operator- done" << std::endl; )
  return m;
}

#include "Matrix_impl.hpp"

#endif
