#ifndef FARRAY_H
#define FARRAY_H

#include <cstddef>
#include <initializer_list>
 
template <typename T>
class FArray1D {

  private:
    bool ownsData;
    int start_i;

  public:
    size_t ni;
    T* data;

    // Constructor to allocate memory
    FArray1D(size_t ni, int start_i=1)
        : ni(ni), start_i(start_i), data(new T[ni]), ownsData(true) {}

    // Constructor to use existing Fortran array
    FArray1D(T* fortranArray, size_t ni, int start_i=1)
        : ni(ni), start_i(start_i), data(fortranArray), ownsData(false) {}

    // Constructor using an initializer list
    FArray1D(std::initializer_list<T> values, size_t ni, int start_i = 1)
        : ni(ni), start_i(start_i), data(new T[ni]), ownsData(true) {
        size_t idx = 0;
        for (T value : values) {
            if (idx < ni) {
                data[idx++] = value;
            }
        }
    }

    // Destructor to free memory
    ~FArray1D() {
        if (data != nullptr && ownsData) {
            delete[] data;
        }
    }

    // Access operator
    T& operator()(size_t i) {
        return data[i-start_i];
    }

    // Const access operator
    const T& operator()(size_t i) const {
        return data[i-start_i];
    }

    // Method to attach an existing Fortran array
    void attach(T* fortranArray) {
        data = fortranArray;
        ownsData = false;
    }

    void fill(const T& value) {
        for (size_t i=start_i; i <ni+start_i; i++) {
            *this(i) = value;
        }
    }

};


template <typename T>
class FArray2D {

  private:
    bool ownsData;
    int start_i, start_j;

  public:
    size_t ni, nj;
    T* data;

    // Constructor to allocate memory
    FArray2D(size_t ni, size_t nj, int start_i=1, int start_j=1)
        : ni(ni), nj(nj), start_i(start_i), start_j(start_j), data(new T[ni * nj]), ownsData(true) {}

    // Constructor to use existing Fortran array
    FArray2D(T* fortranArray, size_t ni, size_t nj, int start_i=1, int start_j=1)
        : ni(ni), nj(nj), start_i(start_i), start_j(start_j), data(fortranArray), ownsData(false) {}

    // Destructor to free memory
    ~FArray2D() {
        if (data != nullptr && ownsData) {
            delete[] data;
        }
    }

    // Access operator
    T& operator()(size_t i, size_t j) {
        return data[i-start_i + ni*(j-start_j)];
    }

    // Const access operator
    const T& operator()(size_t i, size_t j) const {
        return data[i-start_i + ni*(j-start_j)];
    }

    // Method to attach an existing Fortran array
    void attach(T* fortranArray) {
        data = fortranArray;
        ownsData = false;
    }

    void fill(const T& value) {
        for (size_t j=start_j; j<nj+start_j; j++) {
            for (size_t i=start_i; i<ni+start_i; i++) {
                (*this)(i, j) = value;
            }
        }
    }
};


template <typename T>
class FArray3D {

  private:
    bool ownsData;
    int start_i, start_j, start_k;

  public:
    size_t ni, nj, nk;
    T* data;

    // Constructor to allocate memory
    FArray3D(size_t ni, size_t nj, size_t nk, int start_i=1, int start_j=1, int start_k=1)
        : ni(ni), nj(nj), nk(nk), start_i(start_i), start_j(start_j), start_k(start_k), data(new T[ni * nj * nk]), ownsData(true) {}

    // Constructor to use existing Fortran array
    FArray3D(T* fortranArray, size_t ni, size_t nj, size_t nk, int start_i=1, int start_j=1, int start_k=1)
        : ni(ni), nj(nj), nk(nk), start_i(start_i), start_j(start_j), start_k(start_k), data(fortranArray), ownsData(false) {}

    // Destructor to free memory
    ~FArray3D() {
        if (data != nullptr && ownsData) {
            delete[] data;
        }
    }

    // Access operator
    T& operator()(size_t i, size_t j, size_t k) {
        return data[i-start_i + ni*(j-start_j+nj*(k-start_k))];
    }

    // Const access operator
    const T& operator()(size_t i, size_t j, size_t k) const {
        return data[i-start_i + ni*(j-start_j+nj*(k-start_k))];
    }

    // Method to attach an existing Fortran array
    void attach(T* fortranArray) {
        data = fortranArray;
        ownsData = false;
    }

    // Fill method for 3D array (reverse order: k, j, i)
    void fill(const T& value) {
        for (size_t k=start_k; k<nk+start_k; k++) {
            for (size_t j=start_j; j<nj+start_j; j++) {
                for (size_t i=start_i; i<ni+start_i; i++) {
                    (*this)(i, j, k) = value;
                }
            }
        }
    }
};

template <typename T>
class FArray4D {

  private:
    bool ownsData;
    int start_i, start_j, start_k, start_l;

  public:
    size_t ni, nj, nk, nl;
    T* data;

    // Constructor to allocate memory
    FArray4D(size_t ni, size_t nj, size_t nk, size_t nl, int start_i=1, int start_j=1, int start_k=1, int start_l=1)
        : ni(ni), nj(nj), nk(nk), nl(nl), start_i(start_i), start_j(start_j), start_k(start_k), start_l(start_l), data(new T[ni*nj*nk*nl]), ownsData(true) {}

    // Constructor to use existing Fortran array
    FArray4D(T* fortranArray, size_t ni, size_t nj, size_t nk, int start_i=1, int start_j=1, int start_k=1, int start_l=1)
        : ni(ni), nj(nj), nk(nk), nl(nl), start_i(start_i), start_j(start_j), start_k(start_k), start_l(start_l), data(fortranArray), ownsData(false) {}

    // Destructor to free memory
    ~FArray4D() {
        if (data != nullptr && ownsData) {
            delete[] data;
        }
    }

    // Access operator
    T& operator()(size_t i, size_t j, size_t k, size_t l) {
        return data[i-start_i + ni*(j-start_j+nj*(k-start_k+nk*(l-start_l)))];
    }

    // Const access operator
    const T& operator()(size_t i, size_t j, size_t k, size_t l) const {
        return data[i-start_i + ni*(j-start_j+nj*(k-start_k+nk*(l-start_l)))];
    }

    // Method to attach an existing Fortran array
    void attach(T* fortranArray) {
        data = fortranArray;
        ownsData = false;
    }
};


#endif
