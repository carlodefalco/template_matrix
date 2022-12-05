// Template alias for ColumnVector,
// repeated here to be available outside
// the Matrix class scope
template<class Scalar, std::size_t Rows>
using ColumnVector = Matrix<Scalar, Rows, 1>;
 
template<class S1, class S2, std::size_t Rows, std::size_t Cols>
auto
solve_helper (const Matrix<S1, Rows, Cols>& M, const ColumnVector<S2, Rows>& b) {

  static_assert (Rows == Cols, "Matrix::solve () only implemented for square matrices");

  using S = decltype (M(0,0) * b(0));
  constexpr std::size_t m = Rows;
  constexpr std::size_t n = Cols;
  std::size_t ii{0}, jj{0}, kk{0};
  S f{0};
  S pivot{0.}, maxpivot{0.};
  std::size_t imaxpivot{0};
    
  ColumnVector<S, Rows> x{b};
    
  // Perform LU factorization
  Matrix<S, Rows, Cols> factorized{M};
  std::array<std::size_t, Rows> p{0};

  ii = 0;
  for (auto &pp : p)
      pp = ii++;
  

  for (ii = 0; ii < m-1; ++ii) {
    maxpivot = factorized(p[ii], ii);
    imaxpivot = ii;
    for (kk = ii+1; kk < m; ++kk)
      if (std::abs (factorized (p[kk], ii)) > std::abs (maxpivot)) {
        maxpivot = factorized(p[kk], ii);
        imaxpivot = kk;
      }

    if (imaxpivot != ii)
      std::swap (p[ii], p[imaxpivot]);

    pivot = factorized(p[ii], ii);
    for (jj = ii+1; jj < m; ++jj) {
      factorized(p[jj],ii) /= pivot;

      for (kk = ii+1; kk < n; ++kk)
        factorized (p[jj], kk) += - factorized (p[ii],kk) * factorized (p[jj],ii);

    }
  }

  // Perform Forward Substitution
  for (ii = 0; ii < m; ++ii) {
    f = x(p[ii]);
    for (kk = 0; kk < ii; ++kk)
      f -= factorized (p[ii], kk) * x(p[kk]);
    x(p[ii]) = f;
  }

  // Do Backward Substitution
  for (jj = 1; jj <= m; ++jj) {
    ii = m - jj;
    f = x(p[ii]);
    for (kk = ii+1; kk < n; ++kk)
      f -= factorized (p[ii], kk) * x(p[kk]);
    x(p[ii]) = f / factorized (p[ii], ii);
  }

  return (x);
}

  
extern "C" {

  // forward declaration of dgemm
  void
  dgemm_ (const char *TRANSA, const char *TRANSB, const int *M,
          const int *N, const int *K, const double *ALPHA, const double *A,
          const int *LDA, const double *B, const int *LDB, const double *BETA,
          double *C, const int *LDC);
  
  // forward declaration of sgemm
  void
  sgemm_ (const char *TRANSA, const char *TRANSB, const int *M,
          const int *N, const int *K, const float *ALPHA, const float *A,
          const int *LDA, const float *B, const int *LDB, const float *BETA,
          float *C, const int *LDC);

  // forward declaration of dgesv
  void
  dgesv_ (const int *N, const int *NRHS, double *A,
          const int *LDA, int* IPIV, double *B,
          const int *LDB, int *INFO);

}

// Specialization for double by double multiplication using openblas
template<std::size_t Rows, std::size_t Cols, std::size_t N>
Matrix<double, Rows, Cols> operator* (const Matrix<double, Rows, N>& m1, const Matrix<double, N, Cols>& m2) {
  DOIFVERBOSE ( std::cout << "calling operator*" << std::endl; )
  Matrix<double, Rows, Cols> mtmp{};
  double *p_mtmp = mtmp.get_data ();
  const double *p_m1 = m1.get_data (), *p_m2 = m2.get_data ();

  char TRANS{'N'};
  int m{Rows};
  int n{Cols};
  int k{N};
  double one{1.};
  double zero{1.};
  
  dgemm_ (&TRANS, &TRANS, &m, &n, &k, &one, p_m1,
          &m, p_m2, &k, &zero, p_mtmp, &m);

  DOIFVERBOSE ( std::cout << "operator* done" << std::endl; )
  return mtmp;
}

// Specialization for float by float multiplication using openblas
template<std::size_t Rows, std::size_t Cols, std::size_t N>
auto operator* (const Matrix<float, Rows, N>& m1, const Matrix<float, N, Cols>& m2) {
  DOIFVERBOSE ( std::cout << "calling operator*" << std::endl; )
  using m_type = Matrix<float, Rows, Cols>;
  m_type mtmp{};
  float *p_mtmp = mtmp.get_data ();
  const float *p_m1 = m1.get_data (), *p_m2 = m2.get_data ();

  char TRANS{'N'};
  int m{Rows};
  int n{Cols};
  int k{N};
  float one{1.};
  float zero{1.};
  
  sgemm_ (&TRANS, &TRANS, &m, &n, &k, &one, p_m1,
          &m, p_m2, &k, &zero, p_mtmp, &m);

  DOIFVERBOSE ( std::cout << "operator* done" << std::endl; )
  return mtmp;
}

// Specialization for solve with double RHS and double LHS
// error: enclosing class templates are not explicitly specialized
/*
template<std::size_t Rows, std::size_t Cols>
template<>
auto
Matrix<double, Rows, Cols>::solve (const ColumnVector<double, Rows>& b)
{
  static_assert (Rows == Cols, "Matrix::solve () only implemented for square matrices");
  
  int m{Rows};
  int ione{1};

  auto A = *this;
  auto x = b;
  
  double *p_lhs = A.get_data ();
  double *p_rhs = x.get_data ();

  std::array<int, Rows> ipiv{0};

  int info{0};
  
  dgesv_ (&m, &ione, double *A,
          &m, &(ipiv[0]), double *x,
          &m, &info);

  return x;
  
};
*/

// Lengthy solution ...
// template<std::size_t Rows>
// class Matrix<double, Rows, Rows> { ... };

// Simpler solution ...
// template<class S, std::size_t Rows>
// ColumnVector<S> helper_solve (const Matrix<S, Rows, Rows>&,
//                               const ColumnVector<double, Rows>&);

 /*
template<std::size_t Rows>
auto
solve_helper (const Matrix<double, Rows, Rows>& M, const ColumnVector<double, Rows>& b)
{
  int m{Rows};
  int ione{1};

  auto A = M;
  auto x = b;
  
  double *p_lhs = A.get_data ();
  double *p_rhs = x.get_data ();

  std::array<int, Rows> ipiv{0};

  int info{0};
  
  dgesv_ (&m, &ione, p_lhs, &m, &(ipiv[0]),
          p_rhs, &m, &info);

  std::cout << info << std::endl;
  return x;
  
};
 */
