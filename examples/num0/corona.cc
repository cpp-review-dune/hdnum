// vektoren.cc
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

using namespace std;
using namespace hdnum;

int main ()
{
  using Number = double;

  // Corona Fälle ab 24.2.2020 von https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland
  Vector<Number> U = {16, 18, 21, 26, 53, 66, 117, 150, 188, 240,
                      400, 639, 795, 902, 1139, 1296, 1567, 2369,
                      3062, 3795, 4838, 6012, 7156, 8198, 10999, 13957,
		      16662, 18610, 22672, 27436, 31554};

  cout << "U=" << U << std::endl;

  gnuplot("U.dat",U);
}
