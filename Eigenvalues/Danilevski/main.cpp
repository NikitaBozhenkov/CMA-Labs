#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

#define matrix std::vector<std::vector<double>>

double ComputePolValue(int n, double x, double* k) {
  double s = 1;
  for(int i = n - 1; i >= 0; i--)
    s = s * x + k[i];
  return s;
}

double dihot(int degree, double edgeNegativ, double edgePositiv, double* kf) {
  double x;
  double check;
  for(;;) {
    x = 0.5 * (edgeNegativ + edgePositiv);
    if (fabs(x - edgeNegativ) < 1e-5) {
      check = edgePositiv;
      break;
    } else if (fabs(x - edgePositiv) < 1e-5) {
      check = edgeNegativ;
      break;
    }
    if (ComputePolValue(degree, x, kf) < 0)edgeNegativ = x;
    else edgePositiv = x;
  }

  double pol[degree];
  for(int i = 0; i < degree; ++i) {
    pol[i] = kf[i];
  }
  double pol_first_der[degree-1];
  for(int i = 0; i < degree-1; ++i) {
    pol_first_der[i] = pol[i+1] * (i+1);
  }
  double pol_second_der[degree-2];
  for(int i = 0; i < degree-2; ++i) {
    pol_second_der[i] = pol_first_der[i+1] * (i+1);
  }

  if(ComputePolValue(degree,check,pol) * ComputePolValue(degree-2,check,pol_second_der) < 1e-8) {
  while(ComputePolValue(degree, x, pol) > 1e-2) {
    x = x - ComputePolValue(degree,x,pol) / ComputePolValue(degree-1,x,pol_first_der);
  }
  return check;
  } else return x;
}

void stepUp(int level, double** A, double** B, int* currentRootsCount) {
  double major = 0;
  for(int i = 0; i < level; i++) {
    double s = fabs(A[level][i]);
    if (s > major)major = s;
  }
  major += 1.0;

  currentRootsCount[level] = 0;

  for(int i = 0; i <= currentRootsCount[level - 1]; i++) {
    int signLeft, signRight;
    double edgeLeft, edgeRight;
    double edgeNegativ, edgePositiv;

    if (i == 0)edgeLeft = -major;
    else edgeLeft = B[level - 1][i - 1];

    double rb = ComputePolValue(level, edgeLeft, A[level]);

    if (rb == 0) {
      B[level][currentRootsCount[level]] = edgeLeft;
      currentRootsCount[level]++;
      continue;
    }

    if (rb > 0)signLeft = 1; else signLeft = -1;
    if (i == currentRootsCount[level - 1])edgeRight = major;
    else edgeRight = B[level - 1][i];

    rb = ComputePolValue(level, edgeRight, A[level]);

    if (rb == 0) {
      B[level][currentRootsCount[level]] = edgeRight;
      currentRootsCount[level]++;
      continue;
    }

    if (rb > 0)signRight = 1; else signRight = -1;

    if (signLeft == signRight)continue;
    if (signLeft < 0) {
      edgeNegativ = edgeLeft;
      edgePositiv = edgeRight;
    } else {
      edgeNegativ = edgeRight;
      edgePositiv = edgeLeft;
    }

    B[level][currentRootsCount[level]] = dihot(level, edgeNegativ, edgePositiv, A[level]);
    currentRootsCount[level]++;
  }
}

void polynomRealRoots(int n, double* kf, double* rootsArray, int& rootsCount) {
  auto A = new double* [n + 1];
  auto B = new double* [n + 1];
  int* currentRootsCount = new int[n + 1];

  for(int i = 1; i <= n; i++) {
    A[i] = new double[i];
    B[i] = new double[i];
  }

  for(int i = 0; i < n; i++)A[n][i] = kf[i] / kf[n];
  for(int i1 = n, i = n - 1; i > 0; i1 = i, i--) {
    for(int j1 = i, j = i - 1; j >= 0; j1 = j, j--) {
      A[i][j] = A[i1][j1] * j1 / i1;
    }
  }
  currentRootsCount[1] = 1;
  B[1][0] = -A[1][0];
  for(int i = 2; i <= n; i++)stepUp(i, A, B, currentRootsCount);
  rootsCount = currentRootsCount[n];
  for(int i = 0; i < rootsCount; i++)rootsArray[i] = B[n][i];
  for(int i = 1; i <= n; i++) {
    delete[]A[i];
    delete[]B[i];
  }
  delete[]A;
  delete[]B;
  delete[]currentRootsCount;
}

std::vector<double> MatrixMulVec(const matrix& A, const std::vector<double>& vec) {
  std::vector<double> res(A.size());
  for(int i = 0; i < A.size(); ++i) {
    double sum = 0;
    for(int j = 0; j < A.size(); ++j) {
      sum += A[i][j] * vec[j];
    }
    res[i] = sum;
  }
  return res;
}

void PrintMatrix(const matrix& m) {
  for(int i = 0; i < m.size(); ++i) {
    for(int j = 0; j < m.size(); ++j) {
      std::cout << m[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

int main() {
  int n;
  std::cin >> n;
  matrix A(n), E(n);
  for(int i = 0; i < n; ++i) {
    A[i].resize(n);
    E[i].resize(n);
    E[i][i] = 1;
    for(int j = 0; j < n; ++j) {
      std::cin >> A[i][j];
    }
  }
  std::vector<int> coeffs;
  coeffs.push_back(n);

///// Danilevski
  for(int i = n - 2; i >= 0; --i) {

    // Swap
    if (std::abs(A[i + 1][i]) < 1e-8) {
      bool swaped = false;
      for(int j = i - 1; j >= 0; --j) {
        if (std::abs(A[i + 1][j]) > 1e-8) {
          std::swap(A[i], A[j]);
          for(int t = 0; t < n; ++t) {
            std::swap(A[t][i], A[t][j]);
            std::swap(E[t][i], E[t][j]);
          }
          swaped = true;
          break;
        }
      }
      if (!swaped) {
        coeffs.push_back(i + 1);
        --i;
        if (i == -1) continue;
      }
    }

    // Dividing & Multiplication
    double div_mul_value = A[i + 1][i];
    for(int j = 0; j < n; ++j) {
      A[j][i] /= div_mul_value;
      E[j][i] /= div_mul_value;
      A[i][j] *= div_mul_value;
    }
    // Subtracting
    for(int j = 0; j < n; ++j) {
      if (j != i) {
        double pm_value = A[i + 1][j];
        for(int t = 0; t < n; ++t) {
          A[t][j] -= pm_value * A[t][i];
          E[t][j] -= pm_value * E[t][i];
        }
        for(int t = 0; t < n; ++t) {
          A[i][t] += pm_value * A[j][t];
        }
      }
    }
  }
  PrintMatrix(A);
  coeffs.push_back(0);
  std::reverse(coeffs.begin(), coeffs.end());

  std::vector<std::vector<double>> polynoms(coeffs.size() - 1);

  for(int i = 0; i < coeffs.size() - 1; ++i) {
    for(int j = coeffs[i]; j < coeffs[i + 1]; ++j) {
      polynoms[i].push_back(A[coeffs[i]][j]);
    }
  }

///// SOLVING EQUATION
  std::map<double, int> roots;
  //Computing roots
  for(auto& polynom : polynoms) {
    int end_index = polynom.size() - 1;
    while (fabs(polynom[end_index]) < 1e-8) {
      ++roots[0];
      --end_index;
    }

    int polynom_order = end_index + 1;
    double polynom_coeffs[polynom_order + 1];
    polynom_coeffs[polynom_order] = 1;
    for(int j = 0; j < polynom_order; ++j) {
      polynom_coeffs[j] = -polynom[polynom_order - j - 1];
    }

    double f_roots[polynom_order];
    int roots_number = 0;
    polynomRealRoots(polynom_order, polynom_coeffs, f_roots, roots_number);

    for(int i = 0; i < roots_number; ++i) {
      bool added = false;
      for(const auto& elem : roots) {
        if (fabs(elem.first - f_roots[i]) < 1e-6) {
          ++roots[elem.first];
          added = true;
          break;
        }
      }
      if (!added) {
        ++roots[f_roots[i]];
      }

    }
  }
  for(const auto& elem : roots) {
    std::cout << "Root - " << elem.first << ", d = " << elem.second;
    if (elem.second == 1 || polynoms.size() == 1) {
      std::vector<double> F_eigenvector(n);
      for(int i = 0; i < n; ++i) {
        F_eigenvector[i] = std::pow(elem.first, n - i - 1);
      }
      std::vector<double> A_eigenvector = MatrixMulVec(E, F_eigenvector);
      std::cout << ", eigenvector: [";
      for(const auto& k : A_eigenvector) {
        std::cout << k << " ";
      }
      std::cout << "]";
    }
    std::cout << std::endl;
  }

  return 0;
}