// zahlen.cc
#include <iostream>

int main ()
{
  unsigned int i; // uninitialisierte natuerliche Zahl
  double x(3.14); // initialisierte Fliessommazahl
  float y(1.0);   // einfache Genauigkeit
  short j(3);     // eine 'kleine' Zahl
  std::cout << "(i+x)*(y+j)=" << (i+x)*(y+j) << std::endl;

  return 0;
}
