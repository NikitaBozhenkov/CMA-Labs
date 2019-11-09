#include <iostream>
#include <vector>

#define Matrix std::vector<std::vector<double>>


int main() {
  size_t order;
  std::cin >> order;
  Matrix matrix(order);
  std::vector<double> B = {
      5, 6, 2, 0, 7, -8, -9, 6, 7, -3
  };

  for(int i = 0; i < order; ++i) {
    matrix[i].resize(4);
  }
  for(int i = 0; i < order; ++i) {
    if (i != 0 && i != order-1) {
      std::cin >> matrix[i][0] >> matrix[i][1] >> matrix[i][2];
    } else if (i == 0){
      std::cin >> matrix[0][1] >> matrix[0][2];
    } else {
      std::cin >> matrix[i][0] >> matrix[i][1];
    }
  }

  // hard-code for the first iteration
//  if(std::abs(matrix[0][0]) < std::abs(matrix[1][0])) {
//    std::swap(matrix[0], matrix[1]);
//  }
//  if (matrix[0][0] - 1 > 0.00001) {
//    double div_value = matrix[0][0];
//    for(auto& elem : matrix[0]) {
//      elem /= div_value;
//    }
//  }
//  for(int i = 0; i < 3; ++i) {
//    matrix[1][i] -= matrix[0][i];
//  }

  for(int i = 0; i < order-1; ++i) {
    // findind max
    if (std::abs(matrix[i][1]) < std::abs(matrix[i+1][0])) {
      // good swap
      for(int j = 1; j < 4; ++j) {
        double temp = matrix[i][j];
        matrix[i][j] = matrix[i+1][j-1];
        matrix[i+1][j-1] = temp;
      }
      std::swap(B[i], B[i+1]);
    }

    //making [i][i]-elem = 1
    double div_value = matrix[i][1];
    B[i] /= div_value;
    for(int j = 1; j < 4; ++j) {
      matrix[i][j] /= div_value;
    }

    //Gauss
    double mul_value = matrix[i+1][0];
    for(int j = 0; j < 3; ++j) {
      matrix[i+1][j] -= matrix[i][j+1]*mul_value;
    }
    B[i+1] -= B[i]*mul_value;
  }

  double div_value = matrix[order-1][1];
  matrix[order-1][1] /= div_value;
  B[order-1] /= div_value;

  //Upper Gauss
  for(int i = order-1; i > 1; --i) {
    for(int j = 2; j < 4; ++j) {
      B[i - j + 1] -= B[i]*matrix[i - j + 1][j];
    }
  }

  B[0] -= matrix[0][2];


  for(const auto& elem : B) {
    std::cout << elem << " ";
  }




  return 0;
}