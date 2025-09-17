#pragma once

#include <myUtils.hpp> // TODO: specialize this to only the utils you need (that is why I split them up in lib)

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
  // This constructor basically makes itself a copy of the input arg
  Array(const Array& other)
    : size_(other.size()), data_(other.raw()) {}
  // For literal
  Array(const std::vector<T>& raw)
    : size_(raw.size()), data_(raw) {}
  
  std::vector<T>& raw() {return data_;}
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
    ++size_;
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
  
  // extends if newRowCount > size_; else return and throw no errors
  void extend(size_t newRowCount) {
    if (newRowCount <= size_)
      return;

    // inefficient
    for(int i=size_; i< newRowCount; ++i)
      (*this).push_back(T{});
    size_ = newRowCount;
  }
  
  void print(int eleStrLen = 5) const;
  std::string toString(int eleStrLen = 5) const;
};

} // namespace MyFem
