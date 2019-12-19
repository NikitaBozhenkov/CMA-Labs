#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

#define matrix std::vector<std::vector<double>>

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

double MaxNorm(const std::vector<double>& vec) {
  double max = -INT32_MAX;
  for(const auto& elem : vec) {
    if (std::abs(elem) > max) max = std::abs(elem);
  }
  return max;
}

double ScalarProduct(const std::vector<double>& vec1, const std::vector<double>& vec2) {
  double res = 0;
  for(int i = 0; i < vec1.size(); ++i) {
    res += vec1[i] * vec2[i];
  }
  return res;
}

std::vector<double> ScalarMulVec(const double& value, std::vector<double> vec) {
  std::vector<double> res(vec.size());
  for(int i = 0; i < res.size(); ++i) {
    res[i] = vec[i] * value;
  }
  return res;
}

std::vector<double> ComputeAvMinusLvNorm(const matrix& A, const std::vector<double>& v, const double& lambda) {
  std::vector<double> Av = MatrixMulVec(A, v);
  std::vector<double> lambda_v = ScalarMulVec(lambda, v);
  std::vector<double> temp_vec(v.size());
  for(int i = 0; i < v.size(); ++i) {
    temp_vec[i] = Av[i] - lambda_v[i];
  }
  return temp_vec;
}

int main() {
  int n;
  std::cin >> n;
  matrix A(n);
  for(int i = 0; i < n; ++i) {
    A[i].resize(n);
    for(int j = 0; j < n; ++j) {
      std::cin >> A[i][j];
    }
  }
  std::vector<double> uKMM(n), uKM(n), uK(n), uKP(n), vKMM(n), vKM(n), vK(n), vKP(n),
      eigen_vec_21(n), eigen_vec_22(n);
  std::vector<std::complex<double>> eigen_vec_31(n), eigen_vec_32(n);
  for(int i = 0; i < n; ++i) {
    uK[i] = rand();
  }
  double lambdaK1, lambdaK2;
  std::complex<double> lambdaK31, lambdaK32;
  bool even_iteration = true;

///// HARD CODE OF 3 ITERATIONS
  for(int i = 0; i < 3; ++i) {
    vKP = MatrixMulVec(A, uK);
    uKP = vKP;
    double norm = MaxNorm(uKP);
    for(int j = 0; j < n; ++j) {
      uKP[j] /= norm;
    }
    uKMM = uKM;
    uKM = uK;
    uK = uKP;
    vKMM = vKM;
    vKM = vK;
    vK = vKP;
  }


///// STARTING CYCLE
  while (true) {
/////// FIRST CASE
    vKP = MatrixMulVec(A, uK);
    uKP = vKP;

    //Normalize u^{k+1}
    double norm = MaxNorm(uKP);
    for(int j = 0; j < n; ++j) {
      uKP[j] /= norm;
    }

    //Finding lambda
    lambdaK1 = ScalarProduct(vKP, uK) / ScalarProduct(uK, uK);

    //Checking break condition
    if (MaxNorm(ComputeAvMinusLvNorm(A, uK, lambdaK1)) < 1e-12) {
      std::cout << "One lambda: " << lambdaK1 << "." << std::endl << "Eigenvector: [";
      for(const auto& elem : uK) {
        std::cout << elem << " ";
      }
      std::cout << "]" << std::endl;
      break;
    }

/////// SECOND CASE
    double vK_norm = MaxNorm(vK);
    lambdaK2 = -INT32_MAX;
    for(int i = 0; i < n; ++i) {
      if (std::abs(uKM[i]) >= 1e-8) {
        double temp_sqrt = sqrt(vKP[i] * vK_norm / uKM[i]);
        if (temp_sqrt > lambdaK2) lambdaK2 = temp_sqrt;
      }
    }

    // Check of even_iteration
    if (even_iteration) lambdaK2 = -lambdaK2;

    // Computing eigenvector
    for(int i = 0; i < n; ++i) {
      eigen_vec_21[i] = vK[i] + lambdaK2 * uKM[i];
      eigen_vec_22[i] = vK[i] - lambdaK2 * uKM[i];
    }

    if (MaxNorm(ComputeAvMinusLvNorm(A, eigen_vec_21, lambdaK2)) < 1e-12) {
      std::cout << "Two sign-opposite lambdas: " << lambdaK2 << " and " << -lambdaK2 << "." << std::endl
                << "First eigenvector: [";
      for(const auto& elem : eigen_vec_21) {
        std::cout << elem << " ";
      }
      std::cout << "]" << std::endl << "Second eigenvector: [";
      for(const auto& elem : eigen_vec_22) {
        std::cout << elem << " ";
      }
      std::cout << "]" << std::endl;;
      break;
    }

/////// THIRD CASE
    double vKM_norm = MaxNorm(vKM);

    // Finding max index
    int max_index = 0;
    double max_value = -INT32_MAX;
    for(int i = 0; i < n; ++i)
      if (fabs(vKP[i]) > max_value) {
        max_value = fabs(vKP[i]);
        max_index = i;
      }

    // Computing r & cos
    double r = sqrt(fabs(
        (vKM[max_index] * vKP[max_index] * vK_norm - vK[max_index] * vK[max_index] * vKM_norm)
            /
                (uKMM[max_index] * vK[max_index] - uKM[max_index] * uKM[max_index] * vKM_norm)
    ));
    double cos = (vKP[max_index] * vK_norm + r * r * uKM[max_index]) / (2 * r * vK[max_index]);

    if (fabs(cos) <= 1) {
      double sin = sqrt(1 - cos * cos);
      lambdaK31.real(r * cos);
      lambdaK31.imag(r * sin);
      lambdaK32.real(r * cos);
      lambdaK32.imag(-r * sin);
      for(int i = 0; i < n; ++i) {
        eigen_vec_31[i] = vK[i] - lambdaK32 * uKM[i];
        eigen_vec_32[i] = vK[i] - lambdaK31 * uKM[i];
      }

      // Computing Au1 & Au2
      std::vector<std::complex<double>> Au1(n), Au2(n);
      for(int i = 0; i < A.size(); ++i) {
        std::complex<double> sum1(0, 0), sum2(0, 0);
        for(int j = 0; j < A.size(); ++j) {
          sum1 += A[i][j] * eigen_vec_31[j];
          sum2 += A[i][j] * eigen_vec_32[j];
        }
        Au1[i] = sum1;
        Au2[i] = sum2;
      }

      // Computing Au - Lu
      double max_norm_1 = -INT32_MAX;
      for(int i = 0; i < n; ++i) {
        double temp_norm = std::abs(Au1[i] - lambdaK31 * eigen_vec_31[i]);
        if (temp_norm >= max_norm_1) max_norm_1 = temp_norm;
      }
      if (max_norm_1 < 1e-8) {
        std::cout << "Complex pair: " << "r = " << r << ", cos = " << cos << "." << std::endl
                  << "First eigenvector: [";
        for(const auto& elem : eigen_vec_31) {
          std::cout << elem << " ";
        }
        std::cout << std::endl << "]" << std::endl << "Second eigenvector: [";
        for(const auto& elem : eigen_vec_32) {
          std::cout << elem << " ";
        }
        std::cout << "]" << std::endl;

        break;
      }
    }

/////// CHANGING STATUSES
    uKMM = uKM;
    uKM = uK;
    uK = uKP;
    vKMM = vKM;
    vKM = vK;
    vK = vKP;
    even_iteration = !even_iteration;
  }

  return 0;
}
