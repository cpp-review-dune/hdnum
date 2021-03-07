// hallohdnum.cc
#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include "hdnum.hh"    // hdnum header

int main ()
{
  hdnum::Vector<float> a(10,3.14); // Feld mit 10 init. Elementen
  a[3] = 1.0;              // Zugriff auf Element 3
    hdnum::DenseMatrix<double> A={{1,1,2},{2,-3,0}};

  std::vector<double>v(A.colsize(),0) ;
  qrhousholder1(A, v);

  for(int i=0;i<A.rowsize();i++){
    for(int j=0;j<A.colsize();j++){
   std::cout<<A(i,j)<<"      ,      ";
    }
    std::cout<<"\n";
  }
   for(int i=0;i<A.colsize();i++){std::cout<<v[i]<<"\n";}
  return 0;
}
