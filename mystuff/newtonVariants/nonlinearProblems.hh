#ifndef NONLINEARFUNCTIONS_HH
#define NONLINEARFUNCTIONS_HH

#include "hdnum.hh"    // hdnum header

using namespace hdnum;

template<typename Vec, typename Matrix>
class MatyasProblem: public UnconstrainedProblem<Vec, Matrix>{
  public:
  Vec objective(const Vec& x) const{
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vec result(2);
    result[0] = 2.0 * 0.26 * x[0] - 0.48 * x[1];
    result[1] = 2.0 * 0.26 * x[1] - 0.48 * x[0];

    return result;
  }
};

template<typename Vec, typename Matrix>
class RosenbrockProblem: public UnconstrainedProblem<Vec, Matrix>{
  public:
    Vec objective(const Vec& x) const{
      if(x.size() != 2){
        HDNUM_ERROR("Size of vector not equal 2");
      }
      Vec result(2);
      result[0] = -2.0 * (1.0 -x[0]) - 400.0 * (x[1] - x[0] * x[0]) * x[0];
      result[1] = 200.0 * (x[1] -x[0] * x[0]);

      return result;
  }
};

template<typename Vec, typename Matrix>
class BohacheskyProblem: public UnconstrainedProblem<Vec, Matrix>{
  public:
    Vec objective(const Vec& x) const{
      if(x.size() != 2){
        HDNUM_ERROR("Size of vector not equal 2");
      }
      Vec result(2);
      result[0] = 2.0 * x[0] + 0.3 * 3.0 * M_PI * sin(3.0 * M_PI * x[0]);
      result[1] = 4.0 * x[1] + 0.4 * 4.0 * M_PI * sin(4.0 * M_PI * x[1]);
      
      return result;
  }
};

template<typename Vec, typename Matrix>
class BoothProblem: public UnconstrainedProblem<Vec, Matrix>{
  public:
    Vec objective(const Vec& x) const{
      if(x.size() != 2){
        HDNUM_ERROR("Size of vector not equal 2");
      }
      Vec result(2);
      result[0] = 2.0 * (x[0] + 2.0 * x[1] - 7.0) + 4.0 * (2.0 * x[0] +x[1]- 5.0);
      result[1] = 4.0 * (x[0] + 2.0 * x[1] - 7.0) + 2.0 * (2.0 * x[0] +x[1]- 5.0);
      return result;  
  }
};

template<typename Vec, typename Matrix>
class BraninProblem: public UnconstrainedProblem<Vec, Matrix>{
  public:
    Vec objective(const Vec& x) const{
      if(x.size() != 2){
        HDNUM_ERROR("Size of vector not equal 2");
      }
      Vec result(2);

      double a = 1.0;
      double b = 5.1 / ( 4.0 * M_PI * M_PI );
      double c = 5.0 / M_PI;

      double r = 6.0;
      double s = 10.0;
      double t = 1.0 / (8 * M_PI);

      result[0] = 2.0 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r) * (-2.0 * b * x[0] + c) - s * (1.0 - t) * sin(x[0]);
      result[1] = 2.0 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r);

      return result;
  }
};

template<typename Vec, typename Matrix>
class ComplexProblem: public UnconstrainedProblem<Vec, Matrix>{
  public:
    Vec objective(const Vec& x) const{
      Vec result(1);
      result[0] = x[0]*x[0] + 2.0*x[0] + 3.0;
      return result;
    }
};

template<typename N> 
CVector<N> complexFunction(const CVector<N>& x){
  CVector<N> result(1);
  result[0] = x[0]*x[0] + 2.0*x[0] + 3.0;
  return result;
}


template<typename N>
class ConstrainedProblem1:public ConstrainedMinimizationProblem<N>{
  public:
  Vector<N> objective(const Vector<N>& x) const
  {
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vector<N> result(1);
    result[0] = (x[0] - 1.0 ) * (x[0] - 1.0 ) + (x[1] - 2.5 ) *  (x[1] - 2.5 );

    return result;
  }

  Vector<N> gradient(const Vector<N>& x) const
  {
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vector<N> result(2);
    result[0] = 2.0 * (x[0] - 1.0 );
    result[1] = 2.0 * (x[1] - 2.5 );

    return result;
  }

  DenseMatrix<N> A() const{
    DenseMatrix<double> constraints(2,2);
    constraints[0][0] = -1;
    constraints[0][1] = 2;
    constraints[1][0] = 1;
    constraints[1][1] = 2;

    return constraints;
  }

  Vector<N> lowerbounds() const{
    Vector<double> lowerbound = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};
    return lowerbound;
  }

  Vector<N> upperbounds() const{
    Vector<double> upperbound = {2,6};
    return upperbound;
  }
};

template<typename N>
class ConstrainedProblem2:public ConstrainedMinimizationProblem<N>{
  public:
  Vector<N> objective(const Vector<N>& x) const
  {
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vector<N> result(1);
    result[0] = x[0]*x[0] + x[1]*x[1] - 4.0 * x[0] - 5.0 * x[1] + 2.0;

    return result;
  }

  Vector<N> gradient(const Vector<N>& x) const
  {
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vector<N> result(2);
    result[0] = 2.0 * x[0] - 4.0;
    result[1] = 2.0 * x[1] - 5.0;

    return result;
  }

  DenseMatrix<N> A() const{
    DenseMatrix<double> constraints(2,2);
    constraints[0][0] = -2;
    constraints[0][1] = -1;
    constraints[1][0] = 1;
    constraints[1][1] = 0;

    return constraints;
  }

  Vector<N> lowerbounds() const{
    Vector<double> lowerbound = {-2, 0};
    return lowerbound;
  }

  Vector<N> upperbounds() const{
    Vector<double> upperbound = {std::numeric_limits<int>::max(),std::numeric_limits<int>::max()};
    return upperbound;
  }
};

template<typename N> 
class ConstrainedProblem3:public ConstrainedMinimizationProblem<N>{
  public:
  Vector<N> objective(const Vector<N>& x) const{
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vector<N> result(1);
    result[0] = 2.0 * x[0]*x[0] + x[1]*x[1] - 2.0 * x[0] * x[1] - 4.0 * x[0] - 6.0 * x[1];

    return result;
  }

  Vector<N> gradient(const Vector<N>& x) const
  {
    if(x.size() != 2){
      HDNUM_ERROR("Size of vector not equal 2");
    }
    Vector<N> result(2);
    result[0] = 4.0 * x[0] - 2.0 * x[1] - 4.0;
    result[1] = 2.0 * x[1] - 2.0 * x[0] - 6.0;

    return result;
  }

  DenseMatrix<N> A() const{
    DenseMatrix<double> constraints(4,2);

    constraints[0][0] = -1;
    constraints[0][1] = 0;
    constraints[1][0] = 0;
    constraints[1][1] = -1;
    constraints[2][0] = 1;
    constraints[2][1] = 1;
    constraints[3][0] = -1;
    constraints[3][1] = 2;

    return constraints;
   }

  Vector<N> lowerbounds() const{
    Vector<double> lowerbound = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};
    return lowerbound;
  }

  Vector<N> upperbounds() const{
    Vector<double> upperbound = {0,0,8,10};
    return upperbound;
  }
};

template<typename N>
class MinDistToCircle: public ConstrainedMinimizationProblem<N>{
  public:
  Vector<N> objective(const Vector<N>& x) const{
    if(x.size() != 3){
      HDNUM_ERROR("Size of vector not equal 3");
    }
    Vector<N> result(1);
    double epsilon = 1e-5;
    result[0] = (x[0]*x[0]+x[1]*x[1])*(1-1/sqrt(x[0]*x[0]+x[1]*x[1]+epsilon))*(1-1/sqrt(x[0]*x[0]+x[1]*x[1]+epsilon)) + (x[2]-10)* (x[2]-10);
    return result;
  }

  Vector<N> gradient(const Vector<N>& x) const
  {
    if(x.size() != 3){
      HDNUM_ERROR("Size of vector not equal 3");
    }
    double epsilon = 1e-5;
    double sum_squared = (x[0]*x[0]+x[1]*x[1]);
    double t = (1-1/sqrt(x[0]*x[0]+x[1]*x[1]+epsilon));
    double z = 1/(sqrt(x[0]*x[0]+x[1]*x[1]+epsilon) * (x[0]*x[0]+x[1]*x[1]+epsilon));
    Vector<N> result(3);
    result[0] = 2.0 * x[0] * t * t + sum_squared * t * z * x[0];
    result[1] = 2.0 * x[1] * t * t + sum_squared * t * z * x[1];
    result[2] = 2.0 * (x[2]-10.0);

    return result;
  }

  DenseMatrix<N> A() const{
    DenseMatrix<double> constraints(4,3);
    constraints[0][0] = 1;
    constraints[0][1] = 1;
    constraints[0][2] = 1;
  
    constraints[1][0] = -1;
    constraints[1][1] = 1;
    constraints[1][2] = 1;
  
    constraints[2][0] = 1;
    constraints[2][1] = -1;
    constraints[2][2] = 1;
    constraints[3][0] = -1;
    constraints[3][1] = -1;
    constraints[3][2] = 1;

    return constraints;
   }

  Vector<N> lowerbounds() const{
    Vector<double> lowerbound = {10, 10, 10,10};
    return lowerbound;
  }

  Vector<N> upperbounds() const{
    Vector<double> upperbound = {std::numeric_limits<int>::max(),std::numeric_limits<int>::max(), std::numeric_limits<int>::max(),std::numeric_limits<int>::max()};
    return upperbound;
  }
};


#endif