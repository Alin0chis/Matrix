#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() {
  rows_ = 1;
  cols_ = 1;
  Allocate();
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1) {
    throw std::invalid_argument("Invalid rows or cols arguments");
  }

  rows_ = rows;
  cols_ = cols;
  Allocate();
}

S21Matrix::~S21Matrix() { Deallocate(); }

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  Allocate();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

void S21Matrix::Allocate() {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::Deallocate() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
  }
  delete[] matrix_;
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}

int S21Matrix::GetCols() const { return cols_; }

int S21Matrix::GetRows() const { return rows_; }

void S21Matrix::SetCols(const int cols) {
  if (cols <= 0) {
    throw std::invalid_argument("Invalid cols argument");
  }

  if (cols == cols_) {
    return;
  }

  double** newMatrix = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    newMatrix[i] = new double[cols];
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols && j < cols_; j++) {
      newMatrix[i][j] = matrix_[i][j];
    }
  }

  int rows = rows_;
  Deallocate();
  matrix_ = newMatrix;
  cols_ = cols;
  rows_ = rows;
}

void S21Matrix::SetRows(const int rows) {
  if (rows <= 0) {
    throw std::invalid_argument("Invalid rows argument");
  }

  if (rows == rows_) {
    return;
  }

  double** newMatrix = new double*[rows];
  for (int i = 0; i < rows; i++) {
    newMatrix[i] = new double[cols_];
  }

  for (int i = 0; i < rows && i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      newMatrix[i][j] = matrix_[i][j];
    }
  }

  int cols = cols_;
  Deallocate();
  matrix_ = newMatrix;
  cols_ = cols;
  rows_ = rows;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) noexcept {
  if ((rows_ != other.rows_) || (cols_ != other.cols_)) {
    return false;
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        return false;
      }
    }
  }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if ((rows_ != other.rows_) || (cols_ != other.cols_)) {
    throw std::invalid_argument("SumMatrix: different dimensions");
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if ((rows_ != other.rows_) || (cols_ != other.cols_)) {
    throw std::invalid_argument("SumMatrix: different dimensions");
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument("Invalid sizes of matrices for multiplying");
  }
  S21Matrix res(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      res(i, j) = 0;
      for (int k = 0; k < cols_; k++) {
        res(i, j) += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = res;
}

S21Matrix S21Matrix::Transpose() noexcept {
  S21Matrix new_matrix(cols_, rows_);

  for (int i = 0; i < cols_; ++i) {
    for (int j = 0; j < rows_; ++j) {
      new_matrix.matrix_[i][j] = matrix_[j][i];
    }
  }

  return new_matrix;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::domain_error("CalcComplements: matrix must be squared");
  }

  S21Matrix new_matrix(rows_, cols_);

  if (rows_ == 1) {
    new_matrix(0, 0) = 1;
  } else {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        S21Matrix minor = this->Minor(i, j);
        new_matrix.matrix_[i][j] = minor.Determinant() * pow(-1, i + j);
      }
    }
  }

  return new_matrix;
}

S21Matrix S21Matrix::Minor(const int i, const int j) {
  if ((i < 0) || (i > rows_ - 1)) {
    throw std::out_of_range("Minor: i argument out of range");
  }

  if ((j < 0) || (j > cols_ - 1)) {
    throw std::out_of_range("Minor: j argument out of range");
  }

  S21Matrix minor(rows_ - 1, cols_ - 1);

  for (int m = 0; m < rows_; ++m) {
    for (int n = 0; n < cols_; ++n) {
      if ((m < i) && (n < j)) {
        minor(m, n) = matrix_[m][n];
      } else if ((m > i) && (n < j)) {
        minor(m - 1, n) = matrix_[m][n];
      } else if ((m < i) && (n > j)) {
        minor(m, n - 1) = matrix_[m][n];
      } else if ((m > i) && (n > j)) {
        minor(m - 1, n - 1) = matrix_[m][n];
      }
    }
  }

  return minor;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::domain_error("Determinant: matrix must be squared");
  }

  double det = 0;

  if (rows_ == 1) {
    det += matrix_[0][0];
  } else {
    for (int i = 0; i < cols_; ++i) {
      S21Matrix minor = this->Minor(0, i);
      det += minor.Determinant() * pow(-1, i) * matrix_[0][i];
    }
  }

  return det;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = this->Determinant();

  if (fabs(det) < EPS) {
    throw std::domain_error("InverseMatrix: matrix determinant is zero");
  }

  S21Matrix compl_matrix = this->CalcComplements();

  S21Matrix transp_matrix = compl_matrix.Transpose();

  transp_matrix.MulNumber(1 / det);

  return transp_matrix;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double other) {
  MulNumber(other);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) noexcept {
  return this->EqMatrix(other);
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res(*this);
  res += other;
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res(*this);
  res -= other;
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res(*this);
  res *= other;
  return res;
}

S21Matrix S21Matrix::operator*(const double& num) noexcept {
  S21Matrix res(*this);
  res *= num;
  return res;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;
  }

  Deallocate();

  rows_ = other.rows_;
  cols_ = other.cols_;
  Allocate();

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }

  return *this;
}

S21Matrix operator*(const double& num, const S21Matrix& other) {
  S21Matrix res(other);
  res *= num;
  return res;
}

double S21Matrix::operator()(const int i, const int j) const {
  if ((i < 0) || (i > rows_ - 1)) {
    throw std::out_of_range("i argument out of range");
  }

  if ((j < 0) || (j > cols_ - 1)) {
    throw std::out_of_range("j argument out of range");
  }

  return matrix_[i][j];
}

double& S21Matrix::operator()(const int i, const int j) {
  if ((i < 0) || (i > rows_ - 1)) {
    throw std::out_of_range("i argument out of range");
  }

  if ((j < 0) || (j > cols_ - 1)) {
    throw std::out_of_range("j argument out of range");
  }

  double& value = matrix_[i][j];

  return value;
}