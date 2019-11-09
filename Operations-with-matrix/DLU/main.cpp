#include <iostream>
#include <vector>

#define Matrix std::vector<std::vector<double>>

int main() {
  size_t order;
  std::cin >> order;
  Matrix matrix(order);
  Matrix L(order);
  Matrix E(order);
  std::vector<double> B = {
      6, 14, 5
  };
  for(int i = 0; i < order; ++i) {
    matrix[i].resize(order);
    L[i].resize(order);
    E[i].resize(order);
  }

  for(size_t i = 0; i < order; ++i) {
    E[i][i] = 1;
    for(size_t j = 0; j < order; ++j) {
      std::cin >> matrix[i][j];
      //matrix[i][j] = (rand() % 100) * 0.1;
    }
  }

  std::vector<int> this_permutation(order);

  for(int i = 0; i < order; ++i) {
    this_permutation[i] = i;
  }

  for(int i = 0; i < order; ++i) {
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
    if(i != 0) {
      std::swap(L[i], L[max_index]);
    }

    //filling L[i][i] elem
    if (matrix[i][i] != 0) {
      L[i][i] = matrix[i][i];
    } else {
      std::cout << "TRUBA" << std::endl;
      return 0;
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
        for(int t = j; t < order; ++t) {
          matrix[j][t] -= L[j][i] * matrix[i][t];
        }
      }
    }
  }

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
      B[i] /= L[i][i];
      L[i][i] = 1;
      for(int j = i + 1; j < order; ++j) {
        B[j] -= B[i] * L[j][i];
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

  return 0;
}