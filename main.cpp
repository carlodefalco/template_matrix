#include <complex>
#include <chrono>
#include <iostream>
#include <random>
#include <typeinfo>
#include "Matrix.hpp"

constexpr std::size_t BIGN = 1000;

Matrix<double, BIGN, BIGN>
randn () {
  static std::random_device rd;
  Matrix<double, BIGN, BIGN> M;

  double x0 = (rd.max () + rd.min ()) / 2.0;
  double dx = (rd.max () - rd.min ());
  for (std::size_t i = 0; i < BIGN; ++i) 
    for (std::size_t j = 0; j < BIGN; ++j) 
      M(i, j) = (static_cast<double> (rd()) - x0) / dx;

  return M;
};
  
int main () {

  std::cout << "construct random matrices " << std::endl;
  Matrix<double, BIGN, BIGN> m1{randn()}, m2(randn()), m3 = randn(), m4{randn()};
  std::cout << std::endl;

  std::cout << "construct matrix as sum " << std::endl;
  Matrix<double, BIGN, BIGN> m5{m1 + m2};
  std::cout << std::endl;

  std::cout << "construct matrix as product " << std::endl;
  Matrix<double, BIGN, BIGN> m6{m1 * m2};
  std::cout << std::endl;

  
  std::cout << "construct zero matrix " << std::endl;
  Matrix<double, BIGN, BIGN> res{};
  std::cout << std::endl;
  
  auto start = std::chrono::steady_clock::now();
  std::cout << "compute ((m1 * m2) - (m3 * m4))' " << std::endl;
  res = ((m1 * m2) - (m3 * m4)).transpose ();
  std::cout << std::endl;
  

  /*  
  auto start = std::chrono::steady_clock::now();
  std::cout << "compute ((m1 * m2) - (m3 * m4))' " << std::endl;
  Matrix<double, BIGN, BIGN> res = ((m1 * m2) - (m3 * m4)).transpose ();
  std::cout << std::endl;
  */
  
  auto end = std::chrono::steady_clock::now();
  std::cout << "Elapsed time in microseconds: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count ()
            << " us" << std::endl;

  DOIFPRINTRES (
  std::cout << "m1 = [ ..." << std::endl;
  m1.print ();
  std::cout << "]; m2 = [ ..." << std::endl;
  m2.print ();
  std::cout << "]; m3 = [ ..." << std::endl;
  m3.print ();
  std::cout << "]; m4 = [ ..." << std::endl;
  m4.print ();
  std::cout << "]; res = [ ..." << std::endl;
  res.print ();
  std::cout << "];" << std::endl; )
}
