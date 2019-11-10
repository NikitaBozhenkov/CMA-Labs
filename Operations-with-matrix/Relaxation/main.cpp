#include <iostream>
#include <vector>
#include <iomanip>

_inline double AXB_Norm(int n, const std::vector<double>& xi) {
  double max_elem = -INT32_MAX;
  for(int i = 0; i < 3; ++i) {
    double cur_elem = -1;
    if (i == 0) {
      cur_elem += n * xi[0] + (n - 2) * xi[1] + xi[2];
    } else if (i == 2) {
      cur_elem += n * xi[2] + (n-2) * xi[1] + xi[0];
    } else {
      cur_elem += xi[0] + xi[2] + xi[1] * n;
    }
    if (std::abs(cur_elem) > max_elem) {
      max_elem = std::abs(cur_elem);
    }
  }
  return max_elem;
}

int main() {
  int n;
  double w;
  std::cin >> n;
  std::cin >> w;
  std::vector<double> xi(3, 1.0 / n);
  const std::vector<double> g = xi;

  if (n == 1) {
    std::cout << "1" << std::endl;
    return 0;
  } else if (n == 2) {
    std::cout << "x, 1-x" << std::endl;
    return 0;
  }

  //Relaxation
  int iterations = 0;
  while (AXB_Norm(n, xi) > 1e-10) {
    ++iterations;
    for(int i = 0; i <= 2; ++i) {
      if (i == 0) {
        double second_till_last_sum = xi[1] * (n-2) + xi[2];
        xi[i] = (1 - w) * xi[i] + w * (-second_till_last_sum / n + g[i]);
      } else if (i == 2) {
        double first_till_pre_last_sum = xi[1] * (n-2) + xi[0];
        xi[i] = (1 - w) * xi[i] + w * (-first_till_pre_last_sum / n + g[i]);
      } else {
        xi[i] = (1 - w) * xi[i] + w * (-(xi[0] + xi[2]) / n + g[i]);
      }
    }
  }
  std::cout << iterations << std::endl;

  for(int i = 0; i < 3; ++i) {
    std::cout << std::fixed << std::setprecision(13) << xi[i] << " ";
  }

  return 0;
}