// vektoren.cc
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

int main ()
{
  hdnum::Vector<float> x(10);        // Vektor mit 10 float-Elementen
  hdnum::Vector<float> y(10,3.14);  // 10 Elemente initialisiert
  hdnum::Vector<float> a;            // uninitialisierter float-Vektor
  
  x[5] = 5.0; // Zugriff erfolgt wie bei Arrays
  x = 5.0; // Alle Werte sind 5.0 (geht bei Arrays nicht!)
  
  a = x + y; // Wäre bei Arrays auch nicht so einfach
  float d = x*y; //Skalarprodukt!

  return 0;
}
