// gaussseidel.cc
#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include "hdnum.hh"    // hdnum header

int main ()
{
  /// Anfangsvektoren (1,2,3) (4,5,6) (7,8,9)
  hdnum::Vector<float> v1(3);
  fill(v1, (float) 1.0,(float) 1.0);
  hdnum::Vector<float> v2(3);
  fill(v2, (float) 4.0,(float) 1.0);
  hdnum::Vector<float> v3(3);
  fill(v3, (float) 7.0,(float) 1.0);
  hdnum::Vector<float> w1 = v1;
  hdnum::Vector<float> w2;
  hdnum::Vector<float> w3;

  hdnum::Vector<float> proj = w1;
  proj *= (w1*v2)/(w1*w1);

  w2 = v2 - proj;

  proj = w1;
  proj *= w1*v3/(w1*w1);

  w3 = v3 - proj;

  proj = w2;
  proj *= w2*v3/(w2*w2);

  w3 -= proj;

  v1.precision(2);
  std::cout << v1 << std::endl;
  std::cout << v2 << std::endl;
  std::cout << v3 << std::endl;

  std::cout << w1 << std::endl;
  std::cout << w2 << std::endl;
  std::cout << w3 << std::endl;

  std::cout << w1*w2 << std::endl;
  std::cout << w2*w3 << std::endl;
  std::cout << w1*w3 << std::endl;

  return 0;
}
