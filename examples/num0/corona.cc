// vektoren.cc
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

using namespace std;
using namespace hdnum;

// vereinbare Typen
using Number = double;
using Vec = Vector<Number>;
using Mat = DenseMatrix<Number>;

// interpoliere $u_0$, $\lambda$ aus t_0 und t_i
Vec interpolation (const Vec& N, const Vec& t, int i)
{
  Vec p(2); // Ergebnis hat zwei Komponenten
  p[0] = N[0];
  p[1] = (log(N[i])-log(N[0]))/(t[i]-t[0]);
  return p;
}

Vec linear_fit (const Vec& N, const Vec& t)
{
  // Erstelle die Matrix A
  Mat A(N.size(),2);
  for (int i=0; i<N.size(); ++i) {
    A[i][0] = 1.0; A[i][1] = t[i]-t[0];
  }

  // Erstelle die rechte Seite
  Vec logN(N.size());
  for (int i=0; i<N.size(); ++i) logN[i] = log(N[i]);

  // berechne Ausgleichsgerade
  Mat AT(A.transpose()); // transponieren 
  Mat ATA(2,2); ATA.mm(AT,A); // A^T*A
  Vec b(2);
  AT.mv(b,logN); // b = A^T * logN
  Vec p(2);
  linsolve(ATA,p,b);
  return p;
}

Vec nonlinear_fit (const Vec& N, const Vec& t, const Vec& p0)
{
  auto F = [&] (const Vec& p) // eine Lambdafunktion 
  {
    Vec r(2,0.0);
    for (int i=0; i<N.size(); ++i) {
      r[0] += (p[0]*exp(p[1]*(t[i]-t[0]))-N[i])
	        *exp(p[1]*(t[i]-t[0]));
      r[1] += (p[0]*exp(p[1]*(t[i]-t[0]))-N[i])*(t[i]-t[0])
	        *exp(p[1]*(t[i]-t[0]));
    }
    return r;
  };
  auto nlp = getNonlinearProblem(F,p0);

  Newton newton;                 // Ein Newtonobjekt
  newton.set_maxit(50);          // Setze diverse Parameter
  newton.set_verbosity(2);
  newton.set_reduction(1e-10);
  newton.set_abslimit(1e-100);
  newton.set_linesearchsteps(8);

  Vec p(p0);            // Objekt für die Loesung
  newton.solve(nlp,p);  // Berechne Lösung
  return p;
}

int main ()
{
  using Number = double;
  using Vec = Vector<Number>;

  // Corona Faelle Deutschland ab dem 24.2.2020
  // https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland
  Vec N = {16, 18, 21, 26, 53, 66, 117, 150, 188, 240,
      400, 639, 795, 902, 1139, 1296, 1567, 2369,
      3062, 3795, 4838, 6012, 7156, 8198, 10999, 13957,
      16662, 18610, 22672, 27436, 31554, 36508, 42288,
      48582,52547,57298,61913,67366,73522,79696,85778,
      91714,95391,99225,103228,108202,113525,117658,
	   120479,123016,125098,127584,130450,133830,137439};
  Vec t(N.size());
  fill(t,0.0,1.0); // Zeitpunkt in Tagen
  vector<string> tstring = {"24/02/2020","25/02/2020","26/02/2020",
    "27/02/2020","28/02/2020","29/02/2020","01/03/2020","02/03/2020","03/03/2020",
    "04/03/2020","05/03/2020","06/03/2020","07/03/2020","08/03/2020","09/03/2020",
    "10/03/2020","11/03/2020","12/03/2020","13/03/2020","14/03/2020","15/03/2020",
    "16/03/2020","17/03/2020","18/03/2020","19/03/2020","20/03/2020","21/03/2020",
    "22/03/2020","23/03/2020","24/03/2020","25/03/2020","26/03/2020","27/03/2020",
    "28/03/2020","29/03/2020","30/03/2020","31/03/2020","01/04/2020","02/04/2020",
    "03/04/2020","04/04/2020","05/04/2020","06/04/2020","07/04/2020","08/04/2020",
    "09/04/2020","10/04/2020","11/04/2020","12/04/2020","13/04/2020","14/04/2020",
			    "15/04/2020","16/04/2020","17/04/2020","18/04/2020"
  }; // Datum zur Ausgabe
  gnuplot("N.dat",tstring,N); // Ausgabe der Daten

  // Italien https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Italien
  // Vec N = {39, 79, 159, 229, 326, 400, 650, 888, 1128, 1694, 2036, 2502, 3089, 3858, 4636, 5883, 7375,
  // 	   9172, 10149, 12462, 15113, 17660, 21157, 24747, 27980, 31506, 35713, 41035, 47021,
  //          53578, 59138, 63927, 69176, 74386, 80539, 86498};
  
  // Vorhersage experiment
  int m=N.size()-5;
  auto Nsub = N.sub(0,m);
  auto tsub = t.sub(0,m);
  int bis=N.size();

  // (1) Interpolation aus erstem und m-1-tem
  auto p1 = interpolation(N,t,m-1);
  Vec u1(bis);
  for (int i=0; i<u1.size(); ++i) u1[i] = p1[0]*exp(p1[1]*i);
  gnuplot("fit1.dat",tstring,u1);
  cout << "Anfangswert (1)" << p1[0] << endl;
  cout << "Verdopplungszeit (1) " << log(2.0)/p1[1] << " Tage" << endl;

  // (2) Ausgleichsgerade für Logarithmen
  auto p2 = linear_fit(Nsub,tsub);
  Vec u2(bis);
  for (int i=0; i<u2.size(); ++i) u2[i] = exp(p2[0])*exp(p2[1]*i);
  gnuplot("fit2.dat",tstring,u2);
  cout << "Anfangswert (2)" << exp(p2[0]) << std::endl;
  cout << "Verdopplungszeit (2) " << log(2.0)/p2[1] << " Tage" << std::endl;

  // (3) nichtlineares Problem
  auto p3(p2);
  for (int k=2; k<=m; k++)
    {
      cout << "k=" << k << std::endl;
      Nsub = N.sub(0,k);
      tsub = t.sub(0,k);
      p3 = nonlinear_fit(Nsub,tsub,p3);
      Vec u3(bis);
      for (int i=0; i<u3.size(); ++i) u3[i] = p3[0]*exp(p3[1]*i);
      gnuplot("fit3.dat",tstring,u3);
      cout << "Anfangswert (3)" << p3[0] << std::endl;
      cout << "Verdopplungszeit (3) " << log(2.0)/p3[1] << " Tage" << std::endl;
    }

  // moving window analysis
  int w=4; // window size
  Vec verdopplungszeiten4;
  for (int i=0; i<w-1; i++) verdopplungszeiten4.push_back(0);
  for (int i=0; i<N.size()-w+1; i++)
    {
	p3[0] = N[i];
	p3[1] = log(2.0)/2.0;
	Nsub = N.sub(i,w);
	tsub = t.sub(i,w);
	p3 = nonlinear_fit(Nsub,tsub,p3);
	cout << "i=" << i << " T_2=" << log(2.0)/p3[1] << " Tage" << std::endl;
	verdopplungszeiten4.push_back(log(2.0)/p3[1]);
    }
  gnuplot("verdopplungszeiten4.dat",tstring,verdopplungszeiten4);

  w=6; // window size
  Vec verdopplungszeiten6;
  for (int i=0; i<w-1; i++) verdopplungszeiten6.push_back(0);
  for (int i=0; i<N.size()-w+1; i++)
    {
	p3[0] = N[i];
	p3[1] = log(2.0)/2.0;
	Nsub = N.sub(i,w);
	tsub = t.sub(i,w);
	p3 = nonlinear_fit(Nsub,tsub,p3);
	cout << "i=" << i << " T_2=" << log(2.0)/p3[1] << " Tage" << std::endl;
	verdopplungszeiten6.push_back(log(2.0)/p3[1]);
    }
  gnuplot("verdopplungszeiten6.dat",tstring,verdopplungszeiten6);

  // Verdopplungszeiten naiv
  Vec verdopplungszeiten_naiv(N.size());
  for (int k=0; k<N.size(); k++)
    {
      int i_max=0;
      for (int i=0; i<=k; i++) if (N[i]<=0.5*N[k]) i_max = i;
      verdopplungszeiten_naiv[k] = k-i_max;
    }
  gnuplot("verdopplungszeiten_naiv.dat",tstring,verdopplungszeiten_naiv);
}
