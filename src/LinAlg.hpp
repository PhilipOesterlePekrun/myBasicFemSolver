#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <algorithm>


template<typename T>
class Array {
  using size_t = std::size_t;
 private:
  size_t size_;
  std::vector<T> data_;

 public:
 // Default ctor
  Array()
    : size_(0), data_(0) {}
  // Standard constructor
  Array(size_t size)
    : size_(size), data_(size) {}
  // This constructor basically makes a zero vector with the same dimensions as other
  Array(const Array& other)
    : size_(other.size()), data_(other.raw()) {}
  // For literal
  Array(const std::vector<T>& raw)
    : size_(raw.size()), data_(raw) {}
    
  const std::vector<T>& raw() const {return data_;}

  // Access by (i, j)
  T& operator()(size_t i) {
    return data_[i];
  }
  const T& operator()(size_t i) const {
    return data_[i];
  }

  size_t size() const {return size_;}
  
  // Expose some raw std::vector functions
  void resize(size_t newSize) {
    size_ = newSize;
    data_.resize(newSize);
  }
  ///void resize(size_t newSize, double val) {data_.resize(newSize);} //# I think unecessary
  void push_back(T ele) {
    size_++;
    data_.push_back(ele);
  };
  
  void deleteIndices(const std::vector<size_t>& indicesToDelete) {
    std::vector<size_t> keepIndices;
    for (size_t i = 0; i < size_; ++i) {
      if (std::find(indicesToDelete.begin(), indicesToDelete.end(), i) == indicesToDelete.end())
        keepIndices.push_back(i);
    }

    std::vector<T> newData(keepIndices.size());
    for (size_t iNew = 0; iNew < keepIndices.size(); ++iNew) {
      size_t iOld = keepIndices[iNew];
      newData[iNew] = (*this)(iOld);
    }

    size_ = keepIndices.size();
    data_ = std::move(newData);
  }
  void extend(size_t newRowCount) {
    if (newRowCount < size_)
      throw std::invalid_argument("New row count must be >= current row count");

    data_.resize(newRowCount, 0.0);
    size_ = newRowCount;
  }
  
  inline void print(int eleStrLen = 5) const { //# we make inline for now
  std::cout<<"[";
  for(int i=0; i<size_; ++i) {
    std::string tmpStr = std::to_string((*this)(i));
    std::string tmpStr2;
    for(int s=0; s<eleStrLen; ++s) {
      if(s >= tmpStr.length())
        tmpStr2 += " ";
      else
        tmpStr2 +=tmpStr[s];
    }
    std::cout<<tmpStr2;
    if(i<size_-1)
      std::cout<<" ";
  }
  std::cout<<"]^T\n";
  std::cout<<std::endl;
}
};

namespace LinAlg {

class Vectord {
  using size_t = std::size_t;
 private:
  size_t size_;
  std::vector<double> data_;

 public:
  // Default ctor
  Vectord()
    : size_(0), data_(0) {}
  // Standard constructor
  Vectord(size_t size)
    : size_(size), data_(size) {}
  // This constructor just duplicates other
  Vectord(const Vectord& other)
    : size_(other.size()), data_(other.raw()) {}
  // For literal
  Vectord(const std::vector<double>& raw)
    : size_(raw.size()), data_(raw) {}
    
  const std::vector<double>& raw() const {return data_;}

  // Access by (i, j)
  double& operator()(size_t i) {
    return data_[i];
  }
  const double& operator()(size_t i) const {
    return data_[i];
  }

  size_t size() const {return size_;}
  
  // Expose some raw std::vector functions
  void resize(size_t newSize) {
    size_ = newSize;
    data_.resize(newSize);
  }
  ///void resize(size_t newSize, double val) {data_.resize(newSize);} //# I think unecessary
  void push_back(double ele) {
    size_++;
    data_.push_back(ele);
  };
  
  void deleteIndices(const std::vector<size_t>& indicesToDelete) {
    std::vector<size_t> keepIndices;
    for (size_t i = 0; i < size_; ++i) {
      if (std::find(indicesToDelete.begin(), indicesToDelete.end(), i) == indicesToDelete.end())
        keepIndices.push_back(i);
    }

    std::vector<double> newData(keepIndices.size());
    for (size_t iNew = 0; iNew < keepIndices.size(); ++iNew) {
      size_t iOld = keepIndices[iNew];
      newData[iNew] = (*this)(iOld);
    }

    size_ = keepIndices.size();
    data_ = std::move(newData);
  }
  void extend(size_t newRowCount) {
    if (newRowCount < size_)
      throw std::invalid_argument("New row count must be >= current row count");

    data_.resize(newRowCount, 0.0);
    size_ = newRowCount;
  }

  // Non-essential/utility functions
  void scale(double factor) {
    for (auto& val : data_) val *= factor;
  }
  
  double L2Norm() const {
    double sum = 0.0;
    for (const double& val : data_)
      sum += val * val;
    return std::sqrt(sum);
  }
  double maxAbsNorm() const {
    double maxVal = 0.0;
    for (const double& val : data_)
      maxVal = std::max(maxVal, std::abs(val));
    return maxVal;
  }
  
  inline void print(int eleStrLen = 5) const { //# we make inline for now
  std::cout<<"[";
  for(int i=0; i<size_; ++i) {
    std::string tmpStr = std::to_string((*this)(i));
    std::string tmpStr2;
    for(int s=0; s<eleStrLen; ++s) {
      if(s >= tmpStr.length())
        tmpStr2 += " ";
      else
        tmpStr2 +=tmpStr[s];
    }
    std::cout<<tmpStr2;
    if(i<size_-1)
      std::cout<<" ";
  }
  std::cout<<"]^T\n";
  std::cout<<std::endl;
}
};
  
class Matrix2d {
  using size_t = std::size_t;
 private:
  size_t nRows_, nCols_;
  Vectord data_;

 public:
  // Default ctor
  Matrix2d()
    : nRows_(0), nCols_(0), data_() {}
  // Standard constructor
  Matrix2d(size_t rows, size_t cols)
    : nRows_(rows), nCols_(cols), data_(rows * cols) {}
  // This constructor just duplicates other
  Matrix2d(const Matrix2d& other)
    : nRows_(other.nRows()), nCols_(other.nCols()), data_(other.raw()) {}
  // This constructor makes the matrix from a raw flat std::vector and sets the rows and cols
  Matrix2d(size_t rows, size_t cols, const std::vector<double>& flatVect)
    : nRows_(rows), nCols_(cols), data_(flatVect) {}
    
  const Vectord& raw() const {return data_;}

  // Access by (i, j)
  double& operator()(size_t i, size_t j) {
    return data_(i * nCols_ + j);
  }
  const double& operator()(size_t i, size_t j) const {
    return data_(i * nCols_ + j);
  }

  size_t nRows() const {return nRows_;}
  size_t nCols() const {return nCols_;}
  
  Vectord rowAt(size_t i) {
    Vectord row;
    for(int j=0; j<nCols_; ++j)
      row.push_back((*this)(i, j));
    return row;
  }
  const Vectord rowAt(size_t i) const {
    Vectord row;
    for(int j=0; j<nCols_; ++j)
      row.push_back((*this)(i, j));
    return row;
  }
  Vectord colAt(size_t j) {
    Vectord col;
    for(int i=0; i<nRows_; ++i)
      col.push_back((*this)(i, j));
    return col;
  }
  const Vectord colAt(size_t j) const {
    Vectord col;
    for(int i=0; i<nRows_; ++i)
      col.push_back((*this)(i, j));
    return col;
  }
  
  void deleteRows(const std::vector<size_t>& rowsToDelete) {
    std::vector<size_t> keepRows;
    for (size_t i = 0; i < nRows_; ++i) {
      if (std::find(rowsToDelete.begin(), rowsToDelete.end(), i) == rowsToDelete.end())
        keepRows.push_back(i);
    }

    Vectord newData(keepRows.size() * nCols_);
    for (size_t iNew = 0; iNew < keepRows.size(); ++iNew) {
      size_t iOld = keepRows[iNew];
      for (size_t j = 0; j < nCols_; ++j)
        newData(iNew * nCols_ + j) = (*this)(iOld, j);
    }

    nRows_ = keepRows.size();
    data_ = std::move(newData);
  }
  void deleteCols(const std::vector<size_t>& colsToDelete) {
    std::vector<size_t> keepCols;
    for (size_t j = 0; j < nCols_; ++j) {
      if (std::find(colsToDelete.begin(), colsToDelete.end(), j) == colsToDelete.end())
        keepCols.push_back(j);
    }

    Vectord newData(nRows_ * keepCols.size());
    for (size_t i = 0; i < nRows_; ++i)
      for (size_t jNew = 0; jNew < keepCols.size(); ++jNew)
        newData(i * keepCols.size() + jNew) = (*this)(i, keepCols[jNew]);

    nCols_ = keepCols.size();
    data_ = std::move(newData);
  }
  void extendRows(size_t newRowCount) {
    if (newRowCount < nRows_)
      throw std::invalid_argument("New row count must be >= current row count");

    data_.resize(newRowCount * nCols_);
    nRows_ = newRowCount;
  }
  void extendCols(size_t newColCount) {
    if (newColCount < nCols_)
      throw std::invalid_argument("New col count must be >= current col count");

    Vectord newData(nRows_ * newColCount);
    for (size_t i = 0; i < nRows_; ++i)
      for (size_t j = 0; j < nCols_; ++j)
        newData(i * newColCount + j) = (*this)(i, j);

    nCols_ = newColCount;
    data_ = std::move(newData);
  }

  // Non-essential/utility functions
  void scale(double sc) {
    for(int i=0; i<nRows_*nCols_; ++i)
      data_(i) *=sc;
  }
  
  const double trace() const {
    size_t n = std::min(nRows_, nCols_);
    double tr = 0.0;
    for(int i=0; i<n; ++i)
      tr += (*this)(i, i);
    return tr;
  }
  
  double frobeniusNorm() const {
    double sum = 0.0;
    for(int i=0; i<nRows_*nCols_; ++i)
      sum += data_(i) * data_(i);
    return std::sqrt(sum);
  }
  double maxAbs() const {
    double maxVal = 0.0;
    for(int i=0; i<nRows_*nCols_; ++i)
      maxVal = std::max(maxVal, std::abs(data_(i)));
    return maxVal;
  }
  double oneNorm() const {
    double maxColSum = 0.0;
    for (size_t j = 0; j < nCols_; ++j) {
      double colSum = 0.0;
      for (size_t i = 0; i < nRows_; ++i)
        colSum += std::abs((*this)(i, j));
      maxColSum = std::max(maxColSum, colSum);
    }
    return maxColSum;
  }
  double infNorm() const {
    double maxRowSum = 0.0;
    for (size_t i = 0; i < nRows_; ++i) {
      double rowSum = 0.0;
      for (size_t j = 0; j < nCols_; ++j)
        rowSum += std::abs((*this)(i, j));
      maxRowSum = std::max(maxRowSum, rowSum);
    }
    return maxRowSum;
  }
  
  inline void print(int eleStrLen = 5) const {
  for(int i=0; i<nRows_; ++i) {
    std::cout<<"[";
    for(int j=0; j<nCols_; ++j) {
      std::string tmpStr = std::to_string((*this)(i, j));
      std::string tmpStr2;
      // Inefficient as fuck but whatever
      for(int s=0; s<eleStrLen; ++s) {
        if(s >= tmpStr.length() || (*this)(i, j)==0)
          tmpStr2 += " ";
        else
          tmpStr2 +=tmpStr[s];
      }
      std::cout<<tmpStr2;
      if(j<nCols_-1)
        std::cout<<" ";
    }
    std::cout<<"]\n";
  }
  std::cout<<std::endl;
}
};

class Matrix3d {
  // implement when necessary
};

class Matrix4d {
  using size_t = std::size_t;
 private:
  size_t nI_, nJ_, nK_, nL_;
  Vectord data_;

 public:
  // Default ctor
  Matrix4d()
    : nI_(0), nJ_(0), nK_(0), nL_(0), data_() {}
  // Standard constructor
  Matrix4d(size_t I, size_t J, size_t K, size_t L)
    : nI_(I), nJ_(J), nK_(K), nL_(L), data_(I * J * K * L) {}
  // This constructor just duplicates other
  Matrix4d(const Matrix4d& other)
    : nI_(other.nI()), nJ_(other.nJ()), nK_(other.nK()), nL_(other.nL()), data_(other.raw()) {}
  // This constructor makes the matrix from a raw flat std::vector and sets the rows and cols
  Matrix4d(size_t I, size_t J, size_t K, size_t L, const std::vector<double>& flatVect)
    : nI_(I), nJ_(J), nK_(K), nL_(L), data_(flatVect) {}
    
  const Vectord& raw() const {return data_;}

  // Access by (i, j)
  double& operator()(size_t i, size_t j, size_t k, size_t l) {
    return data_(i*nJ_ + j*nK_ + k*nL_ + l);
  }
  const double& operator()(size_t i, size_t j, size_t k, size_t l) const {
    return data_(i*nJ_ + j*nK_ + k*nL_ + l);
  }

  size_t nI() const {return nI_;}
  size_t nJ() const {return nJ_;}
  size_t nK() const {return nK_;}
  size_t nL() const {return nL_;}
  
  /*TODO
  Vectord rowAt(size_t i) {
    Vectord row;
    for(int j=0; j<nCols_; ++j)
      row.push_back((*this)(i, j));
    return row;
  }
  const Vectord rowAt(size_t i) const {
    Vectord row;
    for(int j=0; j<nCols_; ++j)
      row.push_back((*this)(i, j));
    return row;
  }
  Vectord colAt(size_t j) {
    Vectord col;
    for(int i=0; i<nRows_; ++i)
      col.push_back((*this)(i, j));
    return col;
  }
  const Vectord colAt(size_t j) const {
    Vectord col;
    for(int i=0; i<nRows_; ++i)
      col.push_back((*this)(i, j));
    return col;
  }
  
  void deleteRows(const std::vector<size_t>& rowsToDelete) {
    std::vector<size_t> keepRows;
    for (size_t i = 0; i < nRows_; ++i) {
      if (std::find(rowsToDelete.begin(), rowsToDelete.end(), i) == rowsToDelete.end())
        keepRows.push_back(i);
    }

    Vectord newData(keepRows.size() * nCols_);
    for (size_t iNew = 0; iNew < keepRows.size(); ++iNew) {
      size_t iOld = keepRows[iNew];
      for (size_t j = 0; j < nCols_; ++j)
        newData(iNew * nCols_ + j) = (*this)(iOld, j);
    }

    nRows_ = keepRows.size();
    data_ = std::move(newData);
  }
  void deleteCols(const std::vector<size_t>& colsToDelete) {
    std::vector<size_t> keepCols;
    for (size_t j = 0; j < nCols_; ++j) {
      if (std::find(colsToDelete.begin(), colsToDelete.end(), j) == colsToDelete.end())
        keepCols.push_back(j);
    }

    Vectord newData(nRows_ * keepCols.size());
    for (size_t i = 0; i < nRows_; ++i)
      for (size_t jNew = 0; jNew < keepCols.size(); ++jNew)
        newData(i * keepCols.size() + jNew) = (*this)(i, keepCols[jNew]);

    nCols_ = keepCols.size();
    data_ = std::move(newData);
  }
  void extendRows(size_t newRowCount) {
    if (newRowCount < nRows_)
      throw std::invalid_argument("New row count must be >= current row count");

    data_.resize(newRowCount * nCols_);
    nRows_ = newRowCount;
  }
  void extendCols(size_t newColCount) {
    if (newColCount < nCols_)
      throw std::invalid_argument("New col count must be >= current col count");

    Vectord newData(nRows_ * newColCount);
    for (size_t i = 0; i < nRows_; ++i)
      for (size_t j = 0; j < nCols_; ++j)
        newData(i * newColCount + j) = (*this)(i, j);

    nCols_ = newColCount;
    data_ = std::move(newData);
  }*/

  // Non-essential/utility functions
  void scale(double sc) {
    for(int i=0; i<nI_*nJ_*nK_*nL_; ++i)
      data_(i) *=sc;
  }
  
  const double trace() const {
    size_t n = std::min(std::min(nI_, nJ_), std::min(nK_, nL_));
    double tr = 0.0;
    for(int i=0; i<n; ++i)
      tr += (*this)(i, i, i, i);
    return tr;
  }
  
  double frobeniusNorm() const {
    double sum = 0.0;
    for(int i=0; i<nI_*nJ_*nK_*nL_; ++i)
      sum += data_(i) * data_(i);
    return std::sqrt(sum);
  }
  double maxAbs() const {
    double maxVal = 0.0;
    for(int i=0; i<nI_*nJ_*nK_*nL_; ++i)
      maxVal = std::max(maxVal, std::abs(data_(i)));
    return maxVal;
  }
  
  /*TODO
  inline void print(int eleStrLen = 5) const {
  for(int i=0; i<nRows_; ++i) {
    std::cout<<"[";
    for(int j=0; j<nCols_; ++j) {
      std::string tmpStr = std::to_string((*this)(i, j));
      std::string tmpStr2;
      // Inefficient as fuck but whatever
      for(int s=0; s<eleStrLen; ++s) {
        if(s >= tmpStr.length() || (*this)(i, j)==0)
          tmpStr2 += " ";
        else
          tmpStr2 +=tmpStr[s];
      }
      std::cout<<tmpStr2;
      if(j<nCols_-1)
        std::cout<<" ";
    }
    std::cout<<"]\n";
  }
  std::cout<<std::endl;
}
*/
};


// Functions not part of any class

// Basic/misc
inline int diracDelta(int i, int j) {
  if(i==j) return 1;
  else return 0;
}

// Operators linalg
inline Vectord scaleVect(double sc, const Vectord& v) {
  Vectord r(v.size());
  for(int i=0; i<v.size(); ++i)
    r(i) += sc*v(i);
  return r;
}
inline Matrix2d scaleMat2(double sc, const Matrix2d& mat) {
  Matrix2d r(mat.nRows(), mat.nCols());
  for(int i=0; i<mat.nRows(); ++i)
    for(int j=0; j<mat.nCols(); ++j)
      r(i, j) = sc*mat(i, j);
  return r;
}
inline Vectord addVects(const Vectord& v1, const Vectord& v2) {
  Vectord r(v1.size());
  for(int i=0; i<v1.size(); ++i)
    r(i) = v1(i)+v2(i);
  return r;
}
inline Matrix2d addMat2s(const Matrix2d& mat1, const Matrix2d& mat2) {
  Matrix2d r(mat1.nRows(), mat1.nCols());
  for(int i=0; i<mat1.nRows(); ++i)
    for(int j=0; j<mat1.nCols(); ++j)
      r(i, j) = mat1(i, j)+mat2(i, j);
  return r;
}
inline double vectDotProduct(const Vectord& v1, const Vectord& v2) {
  double r = 0.0;
  for(int i=0; i<v1.size(); ++i)
    r += v1(i)*v2(i);
  return r;
}
inline Matrix2d vectDyadicProduct(const Vectord& v1, const Vectord& v2) {
  Matrix2d r(v1.size(), v2.size());
  for(int i=0; i<v1.size(); ++i)
    for(int j=0; j<v2.size(); ++j)
      r(i, j) = v1(i)*v2(j);
  return r;
}
inline Vectord mat2TimesVect(const Matrix2d& mat, const Vectord& v) {
  Vectord r(mat.nRows());
  for(int i=0; i<mat.nRows(); ++i)
    r(i) = vectDotProduct(mat.rowAt(i), v);
  return r;
}
inline Matrix2d matTimesMat2(const Matrix2d& mat1, const Matrix2d& mat2) {//unfinished
  int rNr = mat1.nRows();
  int rNc = mat2.nCols();
  int Ni = mat1.nCols(); //==mat2.row()
  Matrix2d r(rNr, rNc);
  for(int i=0; i<rNr; ++i)
    for(int j=0; j<rNc; ++j) {
      r(i, j) = vectDotProduct(mat1.rowAt(i), mat2.colAt(j));
    }
  return r;
}
///inline double mat2DoubleContractionMat2(const Matrixd& mat1, const Matrixd& mat2) {} // TODO

// etc
inline Matrix2d transposeMat2(const Matrix2d& mat) {
  Matrix2d r(mat.nCols(), mat.nRows());
  for(int i=0; i<mat.nCols(); ++i)
    for(int j=0; j<mat.nRows(); ++j)
      r(i, j) = mat(j, i);
  return r;
}
inline Matrix2d makeMatSymmetric(const Matrix2d& mat) {
  return addMat2s(scaleMat2(0.5, mat), scaleMat2(0.5, transposeMat2(mat)));
}

inline double detMat2(const Matrix2d& mat) {
  if(mat.nRows() != mat.nCols()) db::throwAndExit("mat.nRows() != mat.nCols()");
  else if(mat.nRows() == 2) return mat(0, 0)*mat(1, 1) - mat(1, 0)*mat(0, 1);
  else db::throwAndExit("det for n>2 not yet implemented"); // TODO
  return 0.0; // return who even cares
}

// Only for small matrices
inline Matrix2d invertMat2(Matrix2d mat) {
  if(mat.nRows() != mat.nCols())
    db::throwAndExit("mat.nRows() != mat.nCols()");
  else if(mat.nRows() == 1) {
    Matrix2d r(1, 1,
      {1.0/mat(0,0)});
    return r;
  }
  else if(mat.nRows() == 2) {
    Matrix2d r(2, 2,
    {
      mat(1, 1), -1.0*mat(0, 1),
      -1.0*mat(1, 0), mat(0, 0)
    }
    );
    return scaleMat2(1.0/detMat2(mat), r);
  }
  else db::throwAndExit("invert for n>2 not yet implemented");
  return Matrix2d(0, 0); // return empty I guess whatever man
}

} //namespace LinAlg
