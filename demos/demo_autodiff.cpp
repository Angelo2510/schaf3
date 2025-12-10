#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <autodiff.hpp>

using namespace ASC_ode;

template <typename T>
T func1(T x, T y)
{
  // return x * sin(y);
  return 1e6 + y;
}

int plot_data()
{
  const int max_degree = 5;
  const int num_points = 201; // number of sample points in [-1, 1]
  const double x_min = -1.0;
  const double x_max = 1.0;

  std::ofstream out("legendre_data.txt");
  if (!out)
  {
    std::cerr << "Error: could not open output file.\n";
    return 1;
  }

  // Write header
  out << std::setprecision(15);
  out << "x";
  for (int n = 0; n <= max_degree; ++n)
    out << " P" << n;
  out << "\n";

  std::vector<double> P;

  for (int i = 0; i < num_points; ++i)
  {
    double x = x_min + (x_max - x_min) * i / (num_points - 1);

    // Compute Legendre polynomials up to degree max_degree at point x
    ASC_ode::LegendrePolynomials(max_degree, x, P);

    // Write x and all P_n(x)
    out << x;
    for (int n = 0; n <= max_degree; ++n)
    {
      out << " " << P[n];
    }
    out << "\n";
  }

  out.close();
  std::cout << "Data written to legendre_data.txt\n";
  return 0;
}

int main()
{
  double x = 1, y = 2;
  AutoDiff<2> adx = Variable<0>(x);
  AutoDiff<2> ady = Variable<1>(y);

  std::cout << "adx = " << adx << std::endl;
  std::cout << "ady = " << ady << std::endl;

  AutoDiff<2> prod = adx * ady;
  std::cout << "prod = " << prod << std::endl;

  std::cout << "func1(adx, ady) = " << func1(adx, ady) << std::endl;

  double eps = 1e-8;
  std::cout << "numdiff df/dx = " << (func1(x + eps, y) - func1(x - eps, y)) / (2 * eps) << std::endl;
  std::cout << "numdiff df/dy = " << (func1(x, y + eps) - func1(x, y - eps)) / (2 * eps) << std::endl;
  auto r = adx / ady;
  std::cout << "EXP: " << r << std::endl;

  {
    // we can do second derivatives:
    AutoDiff<1, AutoDiff<1>> addx{Variable<0>(2)};
    std::cout << "addx = " << addx << std::endl;
    // func = x*x
    // func' = 2*x
    // func'' = 2
    std::cout << "addx*addx = " << addx * addx << std::endl;

    std::cout << "sin(addx) = " << sin(adx) << std::endl;
  }
  std::vector<AutoDiff<2>> vec{};
  LegendrePolynomials(3, ady, vec);
  std::cout << vec[0] << std::endl;
  plot_data();

  return 0;
}
