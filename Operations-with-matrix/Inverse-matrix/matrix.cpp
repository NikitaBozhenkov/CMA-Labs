#include "matrix.h"

Matrix::Matrix(size_t order) : data_(order), order_(order) {
  for(size_t i = 0; i < order; ++i) {
    data_[i].resize(order);
  }
}

Matrix::Matrix(const std::vector<std::vector<double>>& vector) : data_(vector), order_(vector.size()) {}

size_t Matrix::GetOrder() const {
  return order_;
}

std::vector<std::vector<double>> Matrix::GetMatrix() const {
  return data_;
}

std::istream& operator>>(std::istream& stream, Matrix& matrix) {
  stream >> matrix.order_;

  matrix.data_.resize(matrix.order_);
  for(size_t i = 0; i < matrix.order_; ++i) {
    matrix.data_[i].resize(matrix.order_);
    for(size_t j = 0; j < matrix.order_; ++j) {
      double temp;
      stream >> temp;
      matrix.data_[i][j] = temp;
    }
  }

  return stream;
}

std::ostream& operator<<(std::ostream& stream, Matrix& matrix) {
  for(size_t i = 0; i < matrix.order_; ++i) {
    for(size_t j = 0; j < matrix.order_; ++j) {
      std::cout << matrix.data_[i][j];
    }
    std::cout << std::endl;
  }

  return stream;
}

std::vector<double> Matrix::operator*(const std::vector<double>& vector) {
  if (this->GetOrder() != vector.size()) {
    throw std::invalid_argument("Matrices have different order");
  }
  std::vector<double> ret(GetOrder());
  for(size_t i = 0; i < GetOrder(); ++i) {
    double sum = 0;
    for(size_t j = 0; j < GetOrder(); ++j) {
      sum += GetMatrix()[i][j] * vector[j];
    }
    ret[i] = sum;
  }
  return ret;
}

std::vector<double> operator*(const std::vector<double>& vector, const Matrix& matrix) {
  if (matrix.GetOrder() != vector.size()) {
    throw std::invalid_argument("Matrices have different order");
  }
  std::vector<double> ret(matrix.GetOrder());
  for(size_t i = 0; i < matrix.GetOrder(); ++i) {
    double sum = 0;
    for(size_t j = 0; j < matrix.GetOrder(); ++j) {
      sum += matrix.GetMatrix()[j][i] * vector[j];
    }
    ret[i] = sum;
  }
  return ret;
}

Matrix Matrix::operator+(const Matrix& matrix) {
  if (this->GetOrder() != matrix.GetOrder()) {
    throw std::invalid_argument("Matrices have different order");
  }
  Matrix m(GetOrder());
  for(size_t i = 0; i < GetOrder(); ++i) {
    for(size_t j = 0; j < GetOrder(); ++j) {
      m.data_[i][j] = this->data_[i][j] + matrix.GetMatrix()[i][j];
    }
  }
  return *this;
}

Matrix& Matrix::operator+=(const Matrix& matrix) {
  if (this->GetOrder() != matrix.GetOrder()) {
    throw std::invalid_argument("Matrices have different order");
  }
  for(size_t i = 0; i < GetOrder(); ++i) {
    for(size_t j = 0; j < GetOrder(); ++j) {
      this->data_[i][j] += matrix.GetMatrix()[i][j];
    }
  }
  return *this;
}

Matrix Matrix::operator-(const Matrix& matrix) {
  if (this->GetOrder() != matrix.GetOrder()) {
    throw std::invalid_argument("Matrices have different order");
  }
  Matrix m(GetOrder());
  for(size_t i = 0; i < GetOrder(); ++i) {
    for(size_t j = 0; j < GetOrder(); ++j) {
      m.data_[i][j] = this->data_[i][j] - matrix.GetMatrix()[i][j];
    }
  }
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& matrix) {
  if (this->GetOrder() != matrix.GetOrder()) {
    throw std::invalid_argument("Matrices have different order");
  }
  for(size_t i = 0; i < GetOrder(); ++i) {
    for(size_t j = 0; j < GetOrder(); ++j) {
      this->data_[i][j] -= matrix.GetMatrix()[i][j];
    }
  }
  return *this;
}