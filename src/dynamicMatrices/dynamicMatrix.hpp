#ifndef DYNAMICMATRIX_HPP
#define DYNAMICMATRIX_HPP
//#include <iostream>
#include "baseDynamicMatrix.hpp"

namespace dynamicMatrices {
  template <typename T>
  class dynamicMatrix : public baseDynamicMatrix<T> {
    using Matrix = baseDynamicMatrix<T>;

  public:
    // dynamicMatrix (size_t dim);
    dynamicMatrix (size_t rows, size_t columns) : Matrix (rows, columns) {}
    dynamicMatrix (const dynamicMatrix &other)
      : Matrix (other) {} /*
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
}*/
    dynamicMatrix (dynamicMatrix &&other) noexcept : Matrix (other) {}
    /*
: ROWS_ (other.ROWS_),
COLUMNS_ (other.COLUMNS_),
addresses_ (other.addresses_),
data_ (other.data_) {
other.addresses_ = nullptr;
other.data_ = nullptr;
}*/
    virtual ~dynamicMatrix () = default; // {}
    dynamicMatrix &operator= (dynamicMatrix other) { this->swap (other); }
    /*
        T *operator[] (size_t I) { return addresses_[I]; }
        T *operator[] (size_t I) const { return addresses_[I]; }
        T *data () { return data_; }
        T *data () const { return data_; }
        T **addresses () { return addresses_; }
        T **addresses () const { return addresses_; }*/
  };
} // namespace dynamicMatrices
#endif // DYNAMICMATRIX_HPP
