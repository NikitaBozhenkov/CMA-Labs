#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

#define Matrix std::vector<std::vector<double>>

_inline void LU_Equation(Matrix& matrix, std::vector<double>& B) {
  int order = matrix.size();
  Matrix L(order);

  for(int i = 0; i < order; ++i) {
    L[i].resize(order);
  }

  for(int i = 0; i < order; ++i) {
    //filling L[i][i] elem
    if (matrix[i][i] != 0) {
      L[i][i] = matrix[i][i];
    } else {
      std::cout << "TRUBA" << std::endl;
      return;
    }

    //Making [i][i]-elem = 1
    for(int j = i; j < order; ++j) {
      matrix[i][j] /= L[i][i];
    }

    //filling L matrix and doing Gauss
    for(int j = i + 1; j < order; ++j) {
      if (matrix[j][i] != 0) {
        L[j][i] = matrix[j][i];
        matrix[j][i] = 0;
        for(int t = i; t < order; ++t) {
          matrix[j][t] -= L[j][i] * matrix[i][t];
        }
      }
    }
  }

  //LY_2 = B'
  for(int i = 0; i < order; ++i) {
    if (std::abs(B[i]) > 0.00001) {
      B[i] /= L[i][i];
      for(int j = i + 1; j < order; ++j) {
        B[j] -= B[i] * L[j][i];
      }
    }
  }

  //UX = B''
  for(int i = order - 1; i >= 0; --i) {
    if (std::abs(B[i]) > 0.00001) {
      for(int j = i - 1; j >= 0; --j) {
        B[j] -= B[i] * matrix[j][i];
      }
    }
  }

}

_inline void LDLT_Equation(Matrix& matrix, std::vector<double>& B) {
  int order = matrix.size();
  std::vector<bool> D(matrix.size()); // default - false

  for(int i = 0; i < order; ++i) {
    double temp;

    // vector D filling
    if (matrix[i][i] < 0) {
      D[i] = true;
      temp = -std::sqrt(-matrix[i][i]);
    } else {
      temp = std::sqrt(matrix[i][i]);
    }

    // Gauss
    for(int j = i + 1; j < order; ++j) {
      double mul_value = matrix[i][j] / matrix[i][i];
      for(int t = j; t < order; ++t) {
        matrix[j][t] -= matrix[i][t] * mul_value;
      }
    }

    // row div by sqrt
    for(int j = i; j < order; ++j) {
      matrix[i][j] /= temp;
    }
  }

  //LY_1 = B
  for(int i = 0; i < order; ++i) {
    if (std::abs(B[i]) > 0.00001) {
      B[i] /= matrix[i][i];
      for(int j = i + 1; j < order; ++j) {
        B[j] -= B[i] * matrix[i][j];
      }
    }
  }

  //DY_2 = B'
  for(int i = 0; i < order; ++i) {
    if (D[i]) {
      B[i] = -B[i];
    }
  }

  //L^TX=B''
  for(int i = order - 1; i >= 0; --i) {
    if (std::abs(B[i]) > 0.00001) {
      B[i] /= matrix[i][i];
      for(int j = i - 1; j >= 0; --j) {
        B[j] -= B[i] * matrix[j][i];
      }
    }
  }
}

int main() {
  size_t order;
  std::cin >> order;
  Matrix matrix(order);

  for(int i = 0; i < order; ++i) {
    matrix[i].resize(order);
  }

  for(size_t i = 0; i < order; ++i) {
    for(size_t j = i; j < order; ++j) {
      //std::cin >> matrix[i][j];
      matrix[i][j] = (rand() % 100) * 0.1;
      if (i != j) {
        matrix[j][i] = matrix[i][j];
      }
    }
  }

  std::vector<double> true_X(order);
  for(int i = 0; i < order; ++i) {
    true_X[i] = (rand() % 100) * 0.1;
  }

  std::vector<double> B(order);
  for(int i = 0; i < order; ++i) {
    double sum = 0;
    for(int j = 0; j < order; ++j) {
      sum += matrix[i][j] * true_X[j];
    }
    B[i] = sum;
  }

  clock_t start_time = clock();

  LDLT_Equation(matrix, B);
  //LU_Equation(matrix, B);

  clock_t end_time = clock();

  long double search_time = (long double) (end_time - start_time) / CLOCKS_PER_SEC;
  std::cout << "Time:" << search_time << std::endl;

  for(int i = 0; i < order; ++i) {
    if (std::abs(true_X[i] - B[i]) > 0.00001) {
      std::cout << "BAD";
    }
  }
  return 0;
}