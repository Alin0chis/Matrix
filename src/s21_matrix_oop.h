#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <math.h>

#include <iostream>

#define EPS 10E-7

class S21Matrix {
 private:
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated

  void Allocate();
  void Deallocate();

 public:
  S21Matrix();  // Базовый конструктор, инициализирующий матрицу некоторой
                // заранее заданной размерностью.
  S21Matrix(int rows, int cols);  // Параметризированный конструктор с
                                  // количеством строк и столбцов.
  S21Matrix(const S21Matrix& other);  // Конструктор копирования.
  S21Matrix(S21Matrix&& other);  // Конструктор переноса.
  ~S21Matrix();                  // Деструктор.

  // Accessors и mutators
  int GetRows() const;
  int GetCols() const;
  void SetRows(const int rows);
  void SetCols(const int cols);

  // Операции над матрицами
  bool EqMatrix(const S21Matrix& other) noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() noexcept;
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
  S21Matrix Minor(const int i, const int j);

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double& num) noexcept;
  bool operator==(const S21Matrix& other) noexcept;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator*=(const S21Matrix& other);
  S21Matrix operator*=(const double other);
  double operator()(const int i, const int j) const;
  double& operator()(const int i, const int j);
};

#endif  // S21_MATRIX_OOP_H