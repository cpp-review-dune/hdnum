// ausgabe.cc
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

int main ()
{
  hdnum::Vector<float> x(10,1.0);
  hdnum::DenseMatrix<float> M(10,10,1.0);

  std::cout << x << std::endl; // Ausgabe von 7 Nachkommastellen (default)
  std::cout << M << std::endl;

  x.iwidth(2);                // Stellen in Indexausgabe
  x.width(20);                // Anzahl Stellen gesamt
  x.precision(16);            // Anzahl Nachkommastellen

  M.iwidth(2);                // Stellen in Indexausgabe
  M.width(20);                // Anzahl Stellen gesamt
  M.precision(16);            // Anzahl Nachkommastellen

  std::cout << x << std::endl; // Ausgabe von mehr Nachkommastellen
  std::cout << M << std::endl;

  return 0;
}
