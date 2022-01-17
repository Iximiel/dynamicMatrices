#ifndef BASEDYNAMICMATRIX_HPP
#define BASEDYNAMICMATRIX_HPP
//#include <iostream>
#include <numeric>
#include <vector>
namespace dynamicMatrices {
  template <typename T>
  class dynamicMatrix {
  public:
    // dynamicMatrix (size_t dim);
    dynamicMatrix (size_t rows, size_t columns)
      : ROWS_ (rows),
        COLUMNS_ (columns),
        addresses_ (new T *[rows]),
        data_ (new T[columns * rows]) {
      for (size_t I = 0; I < ROWS_; ++I) {
        addresses_[I] = data_ + I * COLUMNS_;
      }
    }
    dynamicMatrix (const dynamicMatrix &other)
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
    dynamicMatrix (dynamicMatrix &&other) noexcept
      : ROWS_ (other.ROWS_),
        COLUMNS_ (other.COLUMNS_),
        addresses_ (other.addresses_),
        data_ (other.data_) {
      other.addresses_ = nullptr;
      other.data_ = nullptr;
    }
    virtual ~dynamicMatrix () {
      delete[] data_;
      delete[] addresses_;
    }
    dynamicMatrix &operator= (dynamicMatrix other) {
      if (&other != this) {
        swap (other);
      }
      return *this;
    }
    dynamicMatrix &operator*= (const T A) {
      for (size_t I = 0; I < ROWS_ * COLUMNS_; ++I) {
        data_[I] *= A;
      }
      return *this;
    }

    dynamicMatrix operator* (const T A) {
      dynamicMatrix toret (this->ROWS_, this->COLUMNS_);
      for (size_t I = 0; I < this->ROWS_ * this->COLUMNS_; ++I) {
        toret.data_[I] = this->data_[I] * A;
      }
      return toret;
    }
    template <typename TT>
    friend dynamicMatrix<TT> operator* (const TT A, const dynamicMatrix<TT> &);

    virtual void swap (dynamicMatrix &other) {
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
  dynamicMatrix<T>
  matMul (const dynamicMatrix<T> &A, const dynamicMatrix<T> &B) {
    if (A.Columns () != B.Rows ()) {
      throw "matMul: cannot multiply: the matrices are not mXn nXp";
    }
    dynamicMatrix<T> toReturn (A.Rows (), B.Columns ());
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
  dynamicMatrix<T> Transpose (const dynamicMatrix<T> &A) {
    dynamicMatrix<T> toReturn (A.Columns (), A.Rows ());
    for (size_t I = 0; I < toReturn.Rows (); ++I) {
      for (size_t J = 0; J < toReturn.Columns (); ++J) {
        toReturn[I][J] = A[J][I];
      }
    }
    return toReturn;
  }
  template <typename T>
  dynamicMatrix<T> operator* (const T A, const dynamicMatrix<T> &orig) {
    dynamicMatrix<T> toret (orig.ROWS_, orig.COLUMNS_);
    for (size_t I = 0; I < orig.ROWS_ * orig.COLUMNS_; ++I) {
      toret.data_[I] = orig.data_[I] * A;
    }
    return toret;
  }

  template <typename T>
  dynamicMatrix<T>
  matMulforT (const dynamicMatrix<T> &A, const dynamicMatrix<T> &BTransposed) {
    /// This special multiplication is a row*column multiplication with the
    /// transpose os B, this function takes advantage of the memory alignment to
    /// boost the speed ot the multiplication
    if (A.Columns () != BTransposed.Columns ()) {
      throw "matMul: cannot multiply: the matrices are not mXn nXp";
    }
    dynamicMatrix<T> toReturn (A.Rows (), BTransposed.Rows ());
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
#endif // BASEDYNAMICMATRIX_HPP
