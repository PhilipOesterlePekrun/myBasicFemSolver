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
    : size_(0), data_(0) {};
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
    : size_(0), data_(0) {};
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
  
class Matrixd {
  using size_t = std::size_t;
 private:
  size_t nRows_, nCols_;
  Vectord data_;

 public:
  // Default ctor
  Matrixd()
    : nRows_(0), nCols_(0), data_() {};
  // Standard constructor
  Matrixd(size_t rows, size_t cols)
    : nRows_(rows), nCols_(cols), data_(rows * cols) {}
  // This constructor just duplicates other
  Matrixd(const Matrixd& other)
    : nRows_(other.nRows()), nCols_(other.nCols()), data_(other.raw()) {}
  // This constructor makes the matrix from a raw flat std::vector and sets the rows and cols
  Matrixd(size_t rows, size_t cols, const std::vector<double>& flatVect)
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
inline Matrixd scaleMat(double sc, const Matrixd& mat) {
  Matrixd r(mat.nRows(), mat.nCols());
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
inline Matrixd addMats(const Matrixd& mat1, const Matrixd& mat2) {
  Matrixd r(mat1.nRows(), mat1.nCols());
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
inline Matrixd vectDyadicProduct(const Vectord& v1, const Vectord& v2) {
  Matrixd r(v1.size(), v2.size());
  for(int i=0; i<v1.size(); ++i)
    for(int j=0; j<v2.size(); ++j)
      r(i, j) = v1(i)*v2(j);
  return r;
}
inline Vectord matTimesVect(const Matrixd& mat, const Vectord& v) {
  Vectord r(mat.nRows());
  for(int i=0; i<mat.nRows(); ++i)
    r(i) = vectDotProduct(mat.rowAt(i), v);
  return r;
}
inline Matrixd matTimesMat(const Matrixd& mat1, const Matrixd& mat2) {//unfinished
  int rNr = mat1.nRows();
  int rNc = mat2.nCols();
  int Ni = mat1.nCols(); //==mat2.row()
  Matrixd r(rNr, rNc);
  for(int i=0; i<rNr; ++i)
    for(int j=0; j<rNc; ++j) {
      r(i, j) = vectDotProduct(mat1.rowAt(i), mat2.colAt(j));
    }
  return r;
}
///inline double mat2DoubleContractionMat2(const Matrixd& mat1, const Matrixd& mat2) {} // TODO

// etc
inline Matrixd transposeMat(const Matrixd& mat) {
  Matrixd r(mat.nCols(), mat.nRows());
  for(int i=0; i<mat.nCols(); ++i)
    for(int j=0; j<mat.nRows(); ++j)
      r(i, j) = mat(j, i);
  return r;
}
inline Matrixd makeMatSymmetric(const Matrixd& mat) {
  return addMats(scaleMat(0.5, mat), scaleMat(0.5, transposeMat(mat)));
}

inline double detMat2(const Matrixd& mat) {
  if(mat.nRows() != mat.nCols()) db::throwAndExit("mat.nRows() != mat.nCols()");
  else if(mat.nRows() == 2) return mat(0, 0)*mat(1, 1) - mat(1, 0)*mat(0, 1);
  else db::throwAndExit("det for n>2 not yet implemented"); // TODO
  return 0.0; // return who even cares
}

// Only for small matrices
inline Matrixd invertMat2(Matrixd mat) {
  if(mat.nRows() != mat.nCols()) db::throwAndExit("mat.nRows() != mat.nCols()");
  else if(mat.nRows() == 2) {
    Matrixd r(2, 2,
    {
      mat(1, 1), -1.0*mat(0, 1),
      -1.0*mat(1, 0), mat(0, 0)
    }
    );
    return scaleMat(1.0/detMat2(mat), r);
  }
  else db::throwAndExit("invert for n>2 not yet implemented");
  return Matrixd(0, 0); // return empty I guess whatever man
}

} //namespace LinAlg
