#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>

#define Matrix std::vector<std::vector<double>>

_inline void VectorGaussSub(std::vector<double>& vec1, const std::vector<double>& vec2, int n) {
  if (vec1[n] >= 0.0001) {
    double div = vec1[n] / vec2[n];
    for(int i = 0; i < n; ++i) {
      vec1[i] -= vec2[i] * div;
    }
  }
}

_inline void IdentityGaussSub(std::vector<double>& vec1, const std::vector<double>& vec2, double mul_value) {
  for(int i = 0; i < vec1.size(); ++i) {
    vec1[i] -= vec2[i] * mul_value;
  }
}

_inline void VectorDiv(std::vector<double>& vec1, const double& value) {
  for(auto& elem : vec1) {
    elem /= value;
  }
}

Matrix MatrixMul(const Matrix& a, const Matrix& b) {
  Matrix r(a.size());
  for(auto& i : r) {
    i.resize(a.size());
  }

  for(int i = 0; i < r.size(); ++i) {
    for(int j = 0; j < r.size(); ++j) {
      double elem = 0;
      for(int t = 0; t < r.size(); ++t) {
        elem += a[i][t] * b[t][j];
      }
      r[i][j] = elem;
    }
  }

  return r;
}

int main() {
  size_t order;
  std::cin >> order;
  Matrix matrix(order);
  Matrix E(order);
  for(int i = 0; i < order; ++i) {
    matrix[i].resize(order);
    E[i].resize(order);
  }

  for(size_t i = 0; i < order; ++i) {
    E[i][i] = 1;
    for(size_t j = 0; j < std::min(order, i + 2); ++j) {
      //std::cin >> matrix[i][j];
      matrix[i][j] = (rand() % 100) * 0.1;
    }
    for(size_t j = order - 1; j >= i + 2; --j) {
      matrix[i][j] = 0;
    }
  }
//  Matrix matrix_copy = matrix;

  clock_t start_time = clock();

  for(int i = order - 1; i > 0; --i) {
    //finding max
    if (std::abs(matrix[i][i]) < std::abs(matrix[i - 1][i])) {
      matrix[i].swap(matrix[i - 1]);
      E[i - 1].swap(E[i]);
    }

    // Dividing by [i][i]-elem
    double div_value = matrix[i][i];
    for(int j = i - 1; j >= 0; --j) {
      matrix[i][j] /= div_value;
    }
    for(auto& elem : E[i]) {
      elem /= div_value;
    }

    //Subtracting rows
    double mul_value = matrix[i - 1][i];
    for(int j = i - 1; j >= 0; --j) {
      matrix[i - 1][j] -= matrix[i][j] * mul_value;
    }
    for(int j = 0; j < order; ++j) {
      E[i - 1][j] -= E[i][j] * mul_value;
    }
  }

  for(int i = 0; i < order; ++i) {
    E[0][i] /= matrix[0][0];
  }

  for(int i = 0; i <= order - 1; ++i) {
    if (std::abs(matrix[i][i]) >= 0.00001) {
    } else {
      std::cout << "TRUBA" << std::endl;
      return 0;
    }
    for(int j = i + 1; j < order; ++j) {
      for(int t = 0; t < order; ++t) {
        E[j][t] -= E[i][t] * matrix[j][i];
      }
    }
  }

  clock_t end_time = clock();

  long double search_time = (long double) (end_time - start_time) / CLOCKS_PER_SEC;
  std::cout << search_time << std::endl;

//  Matrix check = MatrixMul(matrix_copy,E);
//  for(int i = 0; i < order; ++i) {
//    if (std::abs(check[i][i] - 1) > 1e-7) {
//      std::cout << "BAD";
//    }
//    for(int j = 0; j < order; ++j) {
//      if (check[i][j] > 1e-7 && j != i) {
//        std::cout << "BAD";
//      }
//    }
//  }
//
//  for(int i = 0; i < order; ++i) {
//    for(int j = 0; j < order; ++j) {
//      std::cout << E[i][j] << " & ";
//    }
//    std::cout << "\\\\" << std::endl;
//  }
  return 0;
}