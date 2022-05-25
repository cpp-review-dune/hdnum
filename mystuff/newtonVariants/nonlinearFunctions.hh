#ifndef NONLINEARFUNCTIONS_HH
#define NONLINEARFUNCTIONS_HH

#include "hdnum.hh"    // hdnum header

using namespace hdnum;

// Functions F(x) which describe the problem F(x) = 0
template<typename N>
Vector<N> functionMatyas(const Vector<N>& x) {
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = 2.0 * 0.26 * x[0] - 0.48 * x[1];
  result[1] = 2.0 * 0.26 * x[1] - 0.48 * x[0];

  return result;
}

template<typename N>
Vector<N> functionRosenbrock(const Vector<N>& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = -2.0 * (1.0 -x[0]) - 400.0 * (x[1] - x[0] * x[0]) * x[0];
  result[1] = 200.0 * (x[1] -x[0] * x[0]);

  return result;
}

template<typename N>
Vector<N> functionBohachesky(const Vector<N>& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = 2.0 * x[0] + 0.3 * 3.0 * M_PI * sin(3.0 * M_PI * x[0]);
  result[1] = 4.0 * x[1] + 0.4 * 4.0 * M_PI * sin(4.0 * M_PI * x[1]);
  
  return result;
}

template<typename N>
Vector<N> functionBooth(const Vector<N>& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = 2.0 * (x[0] + 2.0 * x[1] - 7.0) + 4.0 * (2.0 * x[0] +x[1]- 5.0);
  result[1] = 4.0 * (x[0] + 2.0 * x[1] - 7.0) + 2.0 * (2.0 * x[0] +x[1]- 5.0);
  return result;
}
  
template<typename N>
Vector<N> functionBranin(const Vector<N>& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  using Real = typename std::conditional<std::is_same<std::complex<double>, N>::value, double, N>::type;

  Real a = 1.0;
  Real b = 5.1 / ( 4.0 * M_PI * M_PI );
  Real c = 5.0 / M_PI;

  Real r = 6.0;
  Real s = 10.0;
  Real t = 1.0 / (8 * M_PI);

  result[0] = 2.0 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r) * (-2.0 * b * x[0] + c) - s * (1.0 - t) * sin(x[0]);
  result[1] = 2.0 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r);

  return result;
}

template<typename N>
Vector<N> complexfunc(const Vector<N>& x){
  Vector<N> result(1);
  result[0] = x[0] * x[0] + 2.0*x[0] + 3.0;
  return result;
}

//Functions F(x) which describe the problem min F(x) (min 0.5*||F(x)||^2 for higher dimensional features)
template<typename N> 
Vector<N> functionConstrained1(const Vector<N>& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(1);
  result[0] = (x[0] - 1.0 ) * (x[0] - 1.0 ) + (x[1] - 2.5 ) *  (x[1] - 2.5 );

  return result;
};

template<typename N> 
Vector<N> gradientConstrained1(const Vector<N>& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = 2.0 * (x[0] - 1.0 );
  result[1] = 2.0 * (x[1] - 2.5 );

  return result;
};

template<typename N> 
Vector<N> functionConstrained2(const Vector<N>& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(1);
  result[0] = x[0]*x[0] + x[1]*x[1] - 4.0 * x[0] - 5.0 * x[1] + 2.0;

  return result;
};

template<typename N> 
Vector<N> gradientConstraiend2(const Vector<N>& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = 2.0 * x[0] - 4.0;
  result[1] = 2.0 * x[1] - 5.0;

  return result;
};

template<typename N> 
Vector<N> functionConstrained3(const Vector<N>& x){
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(1);
  result[0] = 2.0 * x[0]*x[0] + x[1]*x[1] - 2.0 * x[0] * x[1] - 4.0 * x[0] - 6.0 * x[1];

  return result;
};

template<typename N> 
Vector<N> gradientConstraiend3(const Vector<N>& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vector<N> result(2);
  result[0] = 4.0 * x[0] - 2.0 * x[1] - 4.0;
  result[1] = 2.0 * x[1] - 2.0 * x[0] - 6.0;

  return result;
};

#endif