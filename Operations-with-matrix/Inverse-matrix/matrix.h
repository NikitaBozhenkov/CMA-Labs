#ifndef INVERSE_MATRIX_MATRIX_H
#define INVERSE_MATRIX_MATRIX_H
#include <iostream>
#include <vector>

class Matrix {
 public:
  explicit Matrix(const std::vector<std::vector<double>>& vector);
  explicit Matrix(size_t order);

  size_t GetOrder() const;
  std::vector<std::vector<double>> GetMatrix() const;

  std::vector<double> operator*(const std::vector<double>& vector);
  friend std::vector<double> operator*(const std::vector<double>& vector, const Matrix& matrix);
  Matrix operator+(const Matrix& matrix);
  Matrix& operator+=(const Matrix& matrix);
  Matrix operator-(const Matrix& matrix);
  Matrix& operator-=(const Matrix& matrix);

  friend std::istream& operator>>(std::istream& stream, Matrix& matrix);
  friend std::ostream& operator<<(std::ostream& stream, Matrix& matrix);

 private:
  std::vector<std::vector<double>> data_;
  size_t order_;
};

#endif //INVERSE_MATRIX_MATRIX_H
