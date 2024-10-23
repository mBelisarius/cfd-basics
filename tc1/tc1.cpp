#include <iostream>

#include "nuenv/core"
#include "nuenv/algorithm"

#include "generator/generator.hpp"
#include "linalg/tridiag.hpp"
#include "numcal/integration.hpp"

using namespace std;
using namespace cfd_basics;

double analyticDiffSolution(double x, double L, double Ta, double Tb,
                            double q_dot_const, double k) {
  return Ta + x * (Tb - Ta) / L + x * (L - x) * q_dot_const / (2.0 * k);
}

double analyticMeanValue(double L, double Ta, double Tb, double q_dot_const, double k) {
  /* Mean value theorem */
  return 0.5 * (Ta + Tb) + q_dot_const * L * L / (12.0 * k);
}

double analyticHeatTransferRate(double x, double L, double Ta, double Tb,
                                double q_dot_const, double k, double A) {
  /* Analytical solution of heat transfer rate for constant cross section area at position `x`
  *  -A * ((-Ta + Tb) * k / L + q_dot_const * (1 - ln(x / L)))
  */
  return (k * (Ta - Tb) / L + q_dot_const * x - q_dot_const * L / 2.0) * A;
}

TridiagMatrix<double> generate(double L, double Ta, double Tb, double q_dot_const, double k, nuenv::Index n) {
  double dx = (L - 0.0) / n;

  nuenv::VectorX<double> xp(n);
  nuenv::VectorX<double> dxp = nuenv::VectorX<double>::Constant(n, dx);
  nuenv::VectorX<double> q_dot(n);
  nuenv::VectorX<double> kp(n);

  nuenv::VectorX<double> dxe(n + 1);
  nuenv::VectorX<double> ke(n + 1);

  xp[0] = dxp[0] / 2;
  for (nuenv::Index i = 1; i < n; i++) {
    xp[i] = xp[i - 1] + (dxp[i - 1] + dxp[i]) / 2;
  }

  for (nuenv::Index i = 0; i < n; i++) {
    q_dot[i] = q_dot_const;
    kp[i] = k;  // kp[i] = k(T[i]);
  }

//  dxe[0] = dxp[0] / 2.0;
//  dxe[n] = dxp[n - 1] / 2.0;
  dxe[0] = dxp[0];
  dxe[n] = dxp[n - 1];
  for (nuenv::Index i = 1; i < n; i++) {
    dxe[i] = xp[i] - xp[i - 1];
  }

  cout << dxp << endl;

  ke[0] = kp[0];  // ke[0] = k(Ta);
  ke[n] = kp[n - 1];  // ke[n] = k(Tb);
  for (nuenv::Index i = 1; i < n; i++) {
    ke[i] = Patankar(kp[i - 1], dxp[i - 1], kp[i], dxe[i]);
  }

  Mesh1D<double> mesh(xp, kp, q_dot);

  auto S = DiscreteSystem<double>();
  S.SetMesh(mesh);
  S.SetBoundaries(Ta, Tb);
  auto M = S.GenerateSystem(dxp, ke, dxe);

  return M;
}

int main() {
  // Construction parameters
  // L, Ta, Tb, q_dot_const, k, A
  double L = 1.0;
  double Ta = 0.0;
  double Tb = 100.0;
  double q_dot_const = 5.0e5;
  double k = 400.0;
  double A = 0.1;

  // Solution parameters
  nuenv::Index N = 10;

  // Variables
  double _dx = (L - 0.0) / N;
  vector<double> dxp(N, _dx);

  vector<double> xp(N);
  xp[0] = dxp[0] / 2;
  for (nuenv::Index i = 1; i < N; i++) {
    xp[i] = xp[i - 1] + (dxp[i - 1] + dxp[i]) / 2;
  }

  vector<double> sol_analytic(N);
  for (nuenv::Index i = 0; i < N; i++) {
    sol_analytic[i] = analyticDiffSolution(xp[i], L, Ta, Tb, q_dot_const, k);
  }

  double T_bar_analytic = analyticMeanValue(L, Ta, Tb, q_dot_const, k);
  double q_analytic0 = analyticHeatTransferRate(0.0, L, Ta, Tb, q_dot_const, k, A);
  double q_analyticL = analyticHeatTransferRate(L, L, Ta, Tb, q_dot_const, k, A);


  auto M = generate(L, Ta, Tb, q_dot_const, k, N);
  auto T = M.Solve();

  vector<double> error(N);
  for (nuenv::Index i = 0; i < N; i++) {
    error[i] = sol_analytic[i] - T[i];
  }

  double T_bar_numeric = 0;
  for (nuenv::Index i = 1; i < N; i++) {
    T_bar_numeric += T[i] * dxp[i] / L;
  }

  double q_numeric0 = -A * k * (T[0] - Ta) / (0.5 * dxp[0]);
  double q_numericL = -A * k * (Tb - T[N - 1]) / (0.5 * dxp[N - 1]);

  std::printf("   p                    xp          Tp_numeric_p         Tp_analytic_p            AbsError_p \n");
  for (nuenv::Index i = 0; i < N; i++) {
    std::printf("%4zu %21.14e %21.14e %21.14e %21.14e \n",
                i, xp[i], T[i], sol_analytic[i], error[i]);
  }

  std::printf(
      "T_bar_analytic = %21.14e \t T_bar_numeric = %21.14e \t error = %21.14e \n",
      T_bar_analytic, T_bar_numeric, T_bar_analytic - T_bar_numeric);
  std::printf(
      "q_analytic(0) = %21.14e \t q_numeric = %21.14e \t error = %21.14e \n",
      q_analytic0, q_numeric0, q_analytic0 - q_numeric0);
  std::printf(
      "q_analytic(L) = %21.14e \t q_numeric = %21.14e \t error = %21.14e \n",
      q_analyticL, q_numericL, q_analyticL - q_numericL);

  return 0;
}
