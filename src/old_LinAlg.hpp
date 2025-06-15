#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <array>


// Static types
template<typename T, std::size_t N>
using array = std::array<T, N>;

template<std::size_t N>
using arrayd = array<double, N>;

template<std::size_t N>
using arrayi = array<int, N>;

template<std::size_t Nr, std::size_t Nc>
using matrixd = array<arrayd<Nc>, Nr>;

// Dynamic types
template<typename T>
using dynArray = std::vector<T>;
using dynArrayi = dynArray<int>;
using dynArrayd = dynArray<double>;
using dynMatrixd = dynArray<dynArrayd>;


// Static functions
namespace LinAlg::Static {
  // Pure linalg functions
  template<std::size_t N>
  inline arrayd<N> scaleVect(double sc, const arrayd<N> v1) {
    arrayd<N> r;
    for(int i=0; i<N; ++i)
      r[i] += sc*v1[i];
    return r;
  }
  template<std::size_t N1r, std::size_t N2c>
  inline matrixd<N1r, N2c> scaleMat(double sc, const matrixd<N1r, N2c> mat) {
    matrixd<N1r, N2c> r;
    for(int i=0; i<N1r; ++i)
      for(int j=0; j<N2c; ++j)
        r[i][j] = sc*mat[i][j];
    return r;
  }
  template<std::size_t N>
  inline double vectDotProduct(const arrayd<N> v1, const arrayd<N> v2) {
    double r = 0;
    for(int i=0; i<N; ++i)
      r += v1[i]*v2[i];
    return r;
  }
  template<std::size_t N1r, std::size_t N2c> // redundant, but N1 is also the num rows so called N1r
  inline matrixd<N1r, N2c> vectDyadicProduct(const arrayd<N1r> v1, const arrayd<N2c> v2) {
    matrixd<N1r, N2c> r;
    for(int i=0; i<N1r; ++i)
      for(int j=0; j<N2c; ++j)
        r[i][j] = v1[i]*v2[j];
    return r;
  }
  
  // Op
  template<std::size_t N>
  inline arrayd<N> vectEleWiseOp(double (*op)(double), const arrayd<N> v) {
    arrayd<N> r;
    for(int i=0; i<N; ++i)
      r[i] = op(v[i]);
    return r;
  }
  template<std::size_t N1r, std::size_t N2c>
  inline matrixd<N1r, N2c> matEleWiseOp(double (*op)(double), const matrixd<N1r, N2c> mat) {
    matrixd<N1r, N2c> r;
    for(int i=0; i<N1r; ++i)
      for(int j=0; i<N2c; ++j)
        r[i][j] += op(mat[i][j]);
    return r;
  }
  
  // Translate between dynamic and static
  template<std::size_t N>
  inline arrayd<N> vectdFromDynVectd(const dynArrayd v) {
    arrayd<N> r;
    for(int i=0; i<N; ++i)
      r[i] = v[i];
    return r;
  }
  template<std::size_t N1r, std::size_t N2c>
  inline matrixd<N1r, N2c> matrixdFromDynMatrixd(const dynMatrixd mat) {
    matrixd<N1r, N2c> r;
    for(int i=0; i<N1r; ++i)
      for(int j=0; j<N2c; ++j)
        r[i][j] = mat[i][j];
    return r;
  }
  
  // Printing
  template<std::size_t N>
  void printVectd(const arrayd<N> v) {
    std::cout<<"[";
    for(int i=0; i<N; ++i)
      std::cout<<v[i]<<" ";
    std::cout<<"]\n";
  }
  template<std::size_t N1r, std::size_t N2c>
  void printMatd(const matrixd<N1r, N2c> mat) {
    for(int i=0; i<N1r; ++i) {
      std::cout<<"[";
      for(int j=0; j<N2c; ++j)
        std::cout<<mat[i][j]<<" ";
      std::cout<<"]\n";
    }
  }
} // namespace LinAlg::Static

// Dynamic functions
namespace LinAlg::Dynamic {
  // Pure linalg functions
  inline dynArrayd scaleVect(double sc, const dynArrayd& v) {
    int size = v.size();
    dynArrayd r(size);
    for(int i=0; i<size; ++i)
      r[i] += sc*v[i];
    return r;
  }
  inline dynMatrixd scaleMat(double sc, const dynMatrixd& mat) {
    dynMatrixd r;
    for(int i=0; i<mat.size(); ++i) {
      r.push_back(dynArrayd(mat[0].size()));
      for(int j=0; j<mat[0].size(); ++j)
        r[i][j] = sc*mat[i][j];
    }
    return r;
  }
  inline dynArrayd vectAdd(const dynArrayd& v1, const dynArrayd& v2) {
    int size = v1.size();
    dynArrayd r(size);
    for(int i=0; i<size; ++i)
      r[i] = v1[i]+v2[i];
    return r;
  }
  inline dynMatrixd matAdd(const dynMatrixd& mat1, const dynMatrixd& mat2) {
    int Nr = mat1.size();
    int Nc = mat1[0].size();
    dynMatrixd r;
    for(int i=0; i<Nr; ++i) {
      r.push_back(dynArrayd(Nc));
      for(int j=0; j<Nc; ++j)
        r[i][j] = mat1[i][j]+mat2[i][j];
    }
    return r;
  }
  inline double vectDotProduct(const dynArrayd& v1, const dynArrayd& v2) {
    int size = v1.size();
    double r;
    for(int i=0; i<size; ++i)
      r += v1[i]*v2[i];
    return r;
  }
  inline dynMatrixd vectDyadicProduct(const dynArrayd& v1, const dynArrayd& v2) {
    dynMatrixd r(v1.size());
    for(int i=0; i<v1.size(); ++i) {
      r.push_back(dynArrayd(v2.size()));
      for(int j=0; j<v2.size(); ++j)
        r[i][j] = v1[i]*v2[j];
    }
    return r;
  }
  inline dynArrayd matTimesVect(const dynMatrixd& mat, const dynArrayd& v) {
    int Nr = mat.size();
    int Nc = mat[0].size();
    dynArrayd r(Nr);
    for(int i=0; i<Nr; ++i)
      r[i] = vectDotProduct(mat[i], v);
    return r;
  }
  inline dynMatrixd matTimesMat(const dynMatrixd& mat1, const dynMatrixd& mat2) {//unfinished
    int rNr = mat1.size();
    int rNc = mat2[0].size();
    int Ni = mat2.size(); //==mat1[0].size()
    dynMatrixd r;
    for(int i=0; i<rNr; ++i) {
      r.push_back(dynArrayd(rNc));
      for(int j=0; j<rNc; ++j) {
        dynArrayd mat2col_j(rNc);
        for(int i2=0; i2<Ni; ++i2)
          mat2col_j[i2] = mat2[i2][j];
        r[i][j] = vectDotProduct(mat1[i], mat2col_j);
      }
    }
    return r;
  }
  
  // Op
  // TODO
  
  // Translate between dynamic and static
  template<std::size_t N>
  inline dynArrayd dynVectdFromVectd(const arrayd<N>& v) {
    dynArrayd r;
    for(int i=0; i<N; ++i)
      r[i] = v[i];
    return r;
  }
  template<std::size_t N1r, std::size_t N2c>
  inline dynMatrixd dynMatrixdFromMatrixd(const matrixd<N1r, N2c>& mat) {
    dynMatrixd r(N1r);
    for(int i=0; i<N1r; ++i) {
      r.push_back(dynArrayd(N2c));
      for(int j=0; j<N2c; ++j)
        r[i][j] = mat[i][j];
    }
    return r;
  }
  
  // Printing
  inline void printVectd(const dynArrayd v) { //# we make inline for now, I dont want to deal with it
    std::cout<<"[";
    for(int i=0; i<v.size(); ++i)
      std::cout<<v[i]<<" ";
    std::cout<<"]\n";
  }
  inline void printMatd(const dynMatrixd mat) {
    for(int i=0; i<mat.size(); ++i) {
      std::cout<<"[";
      for(int j=0; j<mat[0].size(); ++j)
        std::cout<<mat[i][j]<<" ";
      std::cout<<"]\n";
    }
  }
} // namespace LinAlg::Dynamic
