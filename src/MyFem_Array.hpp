#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <algorithm>

namespace MyFem {

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
  
  // Find ele with the value and return the positions of all occurances; if not found, return empty arr
  Array<size_t> find(T val) const {
    Array<size_t> arr = Array<size_t>();
    for(int i=0; i<size_; ++i) {
      if((*this)(i) == val)
        arr.push_back(i);
    }
    return arr;
  }
  
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
  
  void deleteIndices(const Array<size_t>& indicesToDelete) {
    Array<size_t> keepIndices;
    for (size_t i = 0; i < size_; ++i) {
      if(indicesToDelete.find((size_t)i).size()==0)
        keepIndices.push_back(i);
    }

    std::vector<T> newData(keepIndices.size());
    for (size_t iNew = 0; iNew < keepIndices.size(); ++iNew) {
      size_t iOld = keepIndices(iNew);
      newData[iNew] = (*this)(iOld);
    }

    size_ = keepIndices.size();
    data_ = std::move(newData);
  }
  void deleteIndices(const Array<int>& indicesToDelete) { // TODO: think about which we want. size_t or int or both possibly but also maybe not very good to have both
    Array<size_t> keepIndices;
    for (size_t i = 0; i < size_; ++i) {
      if(indicesToDelete.find((size_t)i).size()==0)
        keepIndices.push_back(i);
    }

    std::vector<T> newData(keepIndices.size());
    for (size_t iNew = 0; iNew < keepIndices.size(); ++iNew) {
      size_t iOld = keepIndices(iNew);
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
  
  // TODO: instead of print(), all of these should in general just have a tostring() function so it can also be written to whatever you want, and you can do other operations on the string. print() is unnecessarily limiting
  // TODO: Also this should be a better structure so it looks better
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
  std::cout<<"]^T\n\n";
}
};

} // namespace MyFem
