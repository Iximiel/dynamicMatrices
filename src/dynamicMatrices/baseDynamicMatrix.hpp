#ifndef BASEbaseDynamicMatrix_HPP
#define BASEbaseDynamicMatrix_HPP
//#include <iostream>
#include <numeric>
#include <vector>
namespace dynamicMatrices {
  template <typename T>
  class baseDynamicMatrix {
  public:
    // baseDynamicMatrix (size_t dim);
    baseDynamicMatrix (size_t rows, size_t columns)
      : ROWS_ (rows),
        COLUMNS_ (columns),
        addresses_ (new T *[rows]),
        data_ (new T[columns * rows]) {
      for (size_t I = 0; I < ROWS_; ++I) {
        addresses_[I] = data_ + I * COLUMNS_;
      }
    }
    baseDynamicMatrix (const baseDynamicMatrix &other)
      : ROWS_ (other.ROWS_),
        COLUMNS_ (other.COLUMNS_),
        addresses_ (new T *[other.ROWS_]),
        data_ (new T[other.COLUMNS_ * other.ROWS_]) {
      for (size_t I = 0; I < ROWS_; ++I) {
        addresses_[I] = data_ + I * COLUMNS_;
      }
      for (size_t I = 0; I < ROWS_ * COLUMNS_; ++I) {
        data_[I] = other.data_[I];
      }
    }
    baseDynamicMatrix (baseDynamicMatrix &&other) noexcept
      : ROWS_ (other.ROWS_),
        COLUMNS_ (other.COLUMNS_),
        addresses_ (other.addresses_),
        data_ (other.data_) {
      other.addresses_ = nullptr;
      other.data_ = nullptr;
    }
    virtual ~baseDynamicMatrix () {
      delete[] data_;
      delete[] addresses_;
    }
    baseDynamicMatrix &operator= (baseDynamicMatrix other) {
      this->swap (other);
    }
    virtual void swap (baseDynamicMatrix &other) {
      std::swap (ROWS_, other.ROWS_);
      std::swap (COLUMNS_, other.COLUMNS_);
      std::swap (addresses_, other.addresses_);
      std::swap (data_, other.data_);
    }
    T *operator[] (size_t I) { return addresses_[I]; }
    T *operator[] (size_t I) const { return addresses_[I]; }
    T *data () { return data_; }
    T *data () const { return data_; }
    T **addresses () { return addresses_; }
    T **addresses () const { return addresses_; }
    size_t Rows () const { return ROWS_; }
    size_t Columns () const { return COLUMNS_; }

  private:
    size_t ROWS_{0};
    size_t COLUMNS_{0};
    T **addresses_{nullptr};
    T *data_{nullptr};
  };

  template <typename T>
  baseDynamicMatrix<T>
  matMul (const baseDynamicMatrix<T> &A, const baseDynamicMatrix<T> &B) {
    if (A.Columns () != B.Rows ()) {
      throw "matMul: cannot multiply: the matrices are not mXn nXp";
    }
    baseDynamicMatrix<T> toReturn (A.Rows (), B.Columns ());
    for (size_t I = 0; I < toReturn.Rows (); ++I) {
      for (size_t J = 0; J < toReturn.Columns (); ++J) {
        toReturn[I][J] = T{};
        for (size_t K = 0; K < A.Columns (); ++K) {
          toReturn[I][J] += A[I][K] * B[K][J];
        }
        // std::cerr << toReturn[I][J] << ' ';
      }
      // std::cerr << '\n';
    }
    return toReturn;
  }
  template <typename T>
  baseDynamicMatrix<T> Transpose (const baseDynamicMatrix<T> &A) {
    baseDynamicMatrix<T> toReturn (A.Columns (), A.Rows ());
    for (size_t I = 0; I < toReturn.Rows (); ++I) {
      for (size_t J = 0; J < toReturn.Columns (); ++J) {
        toReturn[I][J] = A[J][I];
      }
    }
    return toReturn;
  }

  template <typename T>
  baseDynamicMatrix<T> matMulforT (
    const baseDynamicMatrix<T> &A, const baseDynamicMatrix<T> &BTransposed) {
    /// This special multiplication is a row*column multiplication with the
    /// transpose os B, this function takes advantage of the memory alignment to
    /// boost the speed ot the multiplication
    if (A.Columns () != BTransposed.Columns ()) {
      throw "matMul: cannot multiply: the matrices are not mXn nXp";
    }
    baseDynamicMatrix<T> toReturn (A.Rows (), BTransposed.Rows ());
    for (size_t I = 0; I < toReturn.Rows (); ++I) {
      for (size_t J = 0; J < toReturn.Columns (); ++J) {
        toReturn[I][J] =
          std::inner_product (A[I], A[I] + A.Columns (), BTransposed[J], T{});
        // std::cerr << toReturn[I][J] << ' ';
      }
      // std::cerr << '\n';
    }
    return toReturn;
  }
} // namespace dynamicMatrices
#endif // BASEbaseDynamicMatrix_HPP
