#include <iostream>
#include <vector>
#include <cmath>

#define Matrix std::vector<std::vector<double>>

int main() {
  size_t order;
  std::cin >> order;
  Matrix matrix(order);
//  Matrix L(order);
  Matrix E(order);
  std::vector<double> B = {
      11, 511, 9841, 87381, 488281, 2015539, 6725601, 19173961, 48427561

  };
  for(int i = 0; i < order; ++i) {
    matrix[i].resize(order);
    E[i].resize(order);
  }

  std::cout << "A" << std::endl;
  for(size_t i = 0; i < order; ++i) {
    E[i][i] = 1;
    for(size_t j = 0; j < order; ++j) {
      matrix[i][j] = std::pow(i+1, j);
      std::cout << matrix[i][j] << " & ";
      //std::cin >> matrix[i][j];
      //matrix[i][j] = (rand() % 100) * 0.1;
    }
    std::cout << "\\\\" << std::endl;
  }

  std::vector<int> this_permutation(order);

  for(int i = 0; i < order; ++i) {
    this_permutation[i] = i;
  }

  for(int i = 0; i < order-1; ++i) {
    //finding max-abs elem
    int max_index = i;
    double max_elem = std::abs(matrix[i][i]);
    for(int j = i + 1; j < order; ++j) {
      double potential_max = std::abs(matrix[j][i]);
      if (potential_max > max_elem) {
        max_elem = potential_max;
        max_index = j;
      }
    }

    // swapping rows
    if (i != max_index) {
      std::swap(matrix[i], matrix[max_index]);
      std::swap(this_permutation[i], this_permutation[max_index]);
    }

    //filling L[i][i] elem
    if (matrix[i][i] == 0) {
      std::cout << "TRUBA" << std::endl;
      return 0;
    }

    //filling L matrix and doing Gauss
    for(int j = i + 1; j < order; ++j) {
      if (matrix[j][i] != 0) {
        double mul_value = matrix[j][i] / matrix[i][i];
        for(int t = i+1; t < order; ++t) {
          matrix[j][t] -= mul_value * matrix[i][t];
        }
      }
    }

    for(int j = i+1; j < order; ++j) {
      matrix[i][j] /= matrix[i][i];
    }
  }

  std::cout << "LU" << std::endl;
  for(int i = 0; i < order; ++i) {
    for(int j = 0; j < order; ++j) {
      std::cout << matrix[i][j] << " & ";
    }
    std::cout << "\\\\" << std::endl;
  }

  std::cout << "Permutations" << std::endl;
  for(int i = 0; i < order; ++i) {
    std::cout << this_permutation[i] << " \\\\ ";
  }
  std::cout << std::endl;

//  DY_1=B
  for(int i = 0; i < order; ++i) {
    if (this_permutation[i] != i) {
      std::swap(B[this_permutation[this_permutation[i]]], B[this_permutation[i]]);
      std::swap(this_permutation[this_permutation[i]], this_permutation[i]);
    }
  }

  //LY_2 = B'
  for(int i = 0; i < order; ++i) {
    if(std::abs(B[i]) > 0.00001) {
      if (matrix[i][i] == 0) {
        std::cout << "Degenerate system";
        return 0;
      }
      B[i] /= matrix[i][i];
      //matrix[i][i] = 1;
      for(int j = i + 1; j < order; ++j) {
        B[j] -= B[i] * matrix[j][i];
      }
    }
  }

  //UX = B''
  for(int i = order-1; i >= 0; --i) {
    if(std::abs(B[i]) > 0.00001) {
      for(int j = i - 1; j >= 0; --j) {
        B[j] -= B[i] * matrix[j][i];
      }
    }
  }

  for(const auto& elem : B) {
    std::cout << elem << " ";
  }

  return 0;
}