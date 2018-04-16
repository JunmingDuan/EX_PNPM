#ifndef __VEC_H_
#define __VEC_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

template<class Type>
class VEC : public std::vector<Type>
{
  private:
    typedef std::vector<Type> Base;
  public:
    VEC() { }

    VEC(int size, Type val = 0) : Base(size, val) { }

    template<class C>
      VEC& operator=(const VEC<C>& vv){
        if((*this).size() != vv.size()){
          std::cerr<<"Different dim"<<std::endl;
          return *this;
        }
        for(size_t i=0;i!=(*this).size();++i){
          (*this)[i]=vv[i];
        }
        return *this;
      }

    VEC& operator+=(const VEC& vv) {
      for(size_t i=0;i!=(*this).size();++i){
        (*this)[i]+=vv[i];
      }
      return *this;
    }

    VEC& operator-=(const VEC& vv) {
      for(size_t i=0;i!=(*this).size();++i){
        (*this)[i]-=vv[i];
      }
      return *this;
    }

    template<class C> VEC& operator*=(const C& CC) {
      for(size_t i=0;i!=(*this).size();++i){
        (*this)[i]*=CC;
      }
      return *this;
    }

    template<class C> VEC& operator/=(const C& CC) {
      for(size_t i=0;i!=(*this).size();++i){
        (*this)[i]/=CC;
      }
      return *this;
    }

    template<class T> friend std::ostream& operator<<(std::ostream&, const VEC&);

    template<class T> friend std::istream& operator>>(std::istream&, VEC&);
};

template<class T>
std::istream& operator>>(std::istream& in, VEC<T>& vv){
  std::cout<<"Please input numbers"<<std::endl;
  T value;
  while(in>>value){
    vv.push_back(value);
  }
  return in;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const VEC<T>& vv){
  //os.precision(4);
  //os<<std::showpos;
  //os.setf(std::ios::scientific);
  for(size_t i=0;i!=vv.size();++i){
    os << vv[i] << " ";
  }

  return os;
}

  template<class T>
VEC<T> operator+(const VEC<T>& v1, const VEC<T>& v2)
{
  if(v1.size() != v2.size()){
    std::cerr<<"diferent dim"<<std::endl;
    return v1;
  }
  VEC<T> v3(v1);
  v3 += v2;
  return v3;
}

  template<class T>
VEC<T> operator-(const VEC<T>& v1, const VEC<T>& v2)
{
  if(v1.size() != v2.size()){
    std::cerr<<"diferent dim"<<std::endl;
    return v1;
  }
  VEC<T> v3(v1);
  v3 -= v2;
  return v3;
}

  template<class T,class T1>
VEC<T> operator*(const VEC<T>& v1, const T1& a)
{
  VEC<T> v3(v1);
  v3 *= a;
  return v3;
}

  template<class T,class T1>
VEC<T> operator*(const T1& a, const VEC<T>& v1)
{
  VEC<T> v3(v1);
  v3 *= a;
  return v3;
}

  template<class T,class T1>
VEC<T> operator/(const VEC<T>& v1, const T1& a)
{
  if(a==0){
    std::cerr<<"Cannot be divided by 0."<<std::endl;
    return v1;
  }
  VEC<T> v3(v1);
  v3 /= a;
  return v3;
}
#endif

