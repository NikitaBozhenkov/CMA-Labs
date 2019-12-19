#include <iostream>
#include <vector>
#include <cmath>

#define matrix std::vector<std::vector<double>>

void PrintMatrix(const matrix& m) {
  for(int i = 0; i < m.size(); ++i) {
    for(int j = 0; j < m.size(); ++j) {
      if (fabs(m[i][j]) < 1e-8) std::cout << 0 << " ";
      else std::cout << m[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

bool CheckStop(const matrix& m) {
  for(int i = 0; i < m.size() - 2; ++i)
    if (fabs(m[i + 1][i]) > 1e-8) {
      double d_sum = m[i][i] * m[i][i] + m[i + 1][i + 1];
      if ((d_sum * d_sum - 4 * (m[i][i] * m[i + 1][i + 1] - m[i][i + 1] * m[i + 1][i])) > 0)
        return false;
    }
  return true;
}

int main() {
  int n;
  std::cin >> n;
  matrix A(n), Q(n);
  for(int i = 0; i < n; ++i) {
    A[i].resize(n);
    Q[i].resize(n);
    Q[i][i] = 1;
    for(int j = 0; j < n; ++j) {
      std::cin >> A[i][j];
    }
  }

  bool symmetric = true;
  // Symmetric check
  for(int i = 0; i < n; ++i)
    for(int j = i + 1; j < n; ++j)
      if (A[i][j] != A[j][i]) symmetric = false;

////// Making Hessenberg
  for(int i = 0; i < n - 2; ++i)
    for(int j = i + 2; j < n; ++j) {
      if (fabs(A[j][i]) < 1e-8) continue;
      // k = i + 1; l = j

      // Computing sin & cos
      double den = sqrt(A[i + 1][i] * A[i + 1][i] + A[j][i] * A[j][i]);
      double cos = A[i + 1][i] / den;
      double sin = A[j][i] / den;

////// QA
      std::vector<double> k_row = A[i + 1];
      std::vector<double> l_row = A[j];
      std::vector<double> Q_k_row = Q[i + 1];
      std::vector<double> Q_l_row = Q[j];

      for(int t = 0; t < n; ++t) {
        A[i + 1][t] = k_row[t] * cos + l_row[t] * sin;
        A[j][t] = -k_row[t] * sin + l_row[t] * cos;
        Q[i + 1][t] = Q_k_row[t] * cos + Q_l_row[t] * sin;
        Q[j][t] = -Q_k_row[t] * sin + Q_l_row[t] * cos;
      }

////// AQ^T
      std::vector<double> col(n);

      for(int t = 0; t < n; t++) {
        col[t] = A[t][i + 1];
        A[t][i + 1] = A[t][i + 1] * cos + A[t][j] * sin;
      }
      for(int t = 0; t < n; t++)
        A[t][j] = -col[t] * sin + A[t][j] * cos;
    }
  PrintMatrix(A);
  std::cout << "Q:" << std::endl;
  PrintMatrix(Q);

////// QR
  while (!CheckStop(A)) {
    std::vector<double> trig_temp;
    for(int i = 0; i < n - 1; ++i) {
      // Computing sin & cos
      double den = sqrt(A[i + 1][i] * A[i + 1][i] + A[i][i] * A[i][i]);
      double cos = A[i][i] / den;
      double sin = -A[i + 1][i] / den;

      // Saving sin & cos
      trig_temp.push_back(cos);
      trig_temp.push_back(sin);

      // Saving A's & Q's rows
      std::vector<double> i_row = A[i];
      std::vector<double> ip_row = A[i + 1];
      std::vector<double> Q_i_row = Q[i];
      std::vector<double> Q_ip_row = Q[i + 1];

      // Rows conversion
      for(int t = 0; t < n; ++t) {
        A[i][t] = i_row[t] * cos - ip_row[t] * sin;
        Q[i][t] = Q_i_row[t] * cos - Q_ip_row[t] * sin;
        A[i + 1][t] = i_row[t] * sin + ip_row[t] * cos;
        Q[i + 1][t] = Q_i_row[t] * sin + Q_ip_row[t] * cos;
      }
    }

    // Reverse multiplication
    for(int i = 0; i < n - 1; ++i) {
      std::vector<double> col(n);
      for(int t = 0; t < n; ++t) {
        col[t] = A[t][i];
        A[t][i] = A[t][i] * trig_temp[2 * i] - A[t][i + 1] * trig_temp[2 * i + 1];
        A[t][i + 1] = col[t] * trig_temp[2 * i + 1] + A[t][i + 1] * trig_temp[2 * i];
      }
    }
  }
  PrintMatrix(A);


///// Finding values and optionally vectors
  for(int i = 0; i < n; ++i) {
    if (i != n - 1 && fabs(A[i + 1][i]) > 1e-8) {
      double a = A[i][i], b = A[i][i + 1], c = A[i + 1][i], d = A[i + 1][i + 1];
      double sqrt_D = sqrt(-((a + d) * (a + d) - 4 * (a * d - b * c)));
      double real = (a + d) / 2;
      double image = sqrt_D / 2;
      std::cout << "Complex pair: (" << real << " +- " << image << "i" << ")" << std::endl;
      ++i;
    } else std::cout << "Eigenvalue: (" << A[i][i] << ")";
    if (symmetric) {
      std::cout << ", eigenvector: [";
      for(int j = 0; j < n; ++j)
        std::cout << Q[i][j] << " ";
      std::cout << "]";
    }
    std::cout << std::endl;
  }

  return 0;
}