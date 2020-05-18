#include <iostream>
#include <math.h>
#include <gmpxx.h>
#include "hdnum.hh"

using namespace hdnum;

template<class T>
void inverse (DenseMatrix<T>& A)
{
  // Vektoren für Äquilibrierung und Permutationen
  int n = A.rowsize();
  Vector<T> s(n);
  Vector<std::size_t> p(n);
  Vector<std::size_t> q(n);

  // LU Zerlegung
  lr_fullpivot(A,p,q);

  // Berechne spaltenweise die Inverse
  Vector<T> x(n);
  Vector<T> b(n);
  DenseMatrix<T> Ainv(n,n);
  for (int i=0; i<n; i++)
    {
      b = T(0.0);
      b[i] = T(1.0);
      permute_forward(p,b);
      solveL(A,b,b);
      solveR(A,x,b);
      permute_backward(q,x);
      Ainv.sc(x,i);
    }
  A = Ainv;
}

template<typename T>
DenseMatrix<double> setupA (const T& integrateK, double alpha, double beta, int n)
{
  DenseMatrix<double> A(n,n);
  double h = (beta-alpha)/n;
  for (int i=0; i<n; i++)
    {
      double ci = alpha+(i+0.5)*h;
      A[i][i] = integrateK(ci,ci,0.5*h);
      for (int j=i+1; j<n; j++)
	A[j][i] = A[i][j] = integrateK(ci,ci+(j-i)*h,0.5*h);
    }
  return A;
}

template<typename T>
Vector<double> setupb (const T& f, double alpha, double beta, int n)
{
  Vector<double> b(n);
  double h = (beta-alpha)/n;
  for (int i=0; i<n; i++)
    b[i] = f(alpha+(i+0.5)*h);
  return b;
}

int main ()
{
  // Die Problemdefinition und die Berechnung der Matrix
  // und rechten Seite erfolgt in double precision

  // das Intervall
  double alpha = 0.0;
  double beta = 1.0;

  // rechte Seite der Integralgleichung
  auto f = [&] (double x) { return 1.0/(1.1+sin(7.0*M_PI*x)); };

  // Integral über den Kern für die Kollokationsmethode
  double gamma = 0.5;
  auto integrateK = [&] (double ci, double cj, double h2) {
    double d = std::abs(ci-cj);
    if (d<h2) return 2.0/(1-gamma)*pow(h2,1-gamma);
    return 1.0/(1-gamma)*(pow(d+h2,1-gamma)-pow(d-h2,1-gamma)); 
  };

  // Diskretisierungsparameter
  int n = 10;

  // Teste Gram-Schmidt Orthogonalisierung
  {
    n = 128;
    gamma = 1e-4;
    DenseMatrix<double> A_double(setupA(integrateK,alpha,beta,n));
    Vector<double> b_double(setupb(f,alpha,beta,n));
    DenseMatrix<double> Q1 = gram_schmidt(A_double);
    DenseMatrix<double> Q1T(Q1.transpose());
    DenseMatrix<double> I1(A_double); I1.mm(Q1T,Q1);
    for (int i=0; i<Q1T.rowsize(); i++) I1[i][i] -= 1.0;
    DenseMatrix<double> Q2 = modified_gram_schmidt(A_double);
    DenseMatrix<double> Q2T(Q2.transpose());
    DenseMatrix<double> I2(A_double); I2.mm(Q2T,Q2);
    for (int i=0; i<Q2T.rowsize(); i++) I2[i][i] -= 1.0;
    std::cout << "cgs " << I1.norm_infty() << " mgs " << I2.norm_infty() << std::endl;
  }
  
  std::vector<int> sizes = {8,32,128,512};
  std::vector<double> gammas = {0.9,0.5,1e-1,1e-2,1e-3,1e-4,1e-5};

  // jetzt rechne mit anderer Genauigkeit
  mpf_set_default_prec(1024);
  //using number = mpf_class;
  using number = float;

  std::vector<number> conditionnumbers;
  std::vector<number> errors;
  
  for (auto gi : gammas)
    for (auto si : sizes)
      {
	n = si;
	gamma = gi;
	
	// Stelle Gleichungssystem auf
	DenseMatrix<double> A_double(setupA(integrateK,alpha,beta,n));
	Vector<double> b_double(setupb(f,alpha,beta,n));

	// kopiere Gleichungssystem
	DenseMatrix<number> A(n,n);
	Vector<number> b(n);
	for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
	    A[i][j] = A_double[i][j];
	for (int i=0; i<n; i++) b[i] = b_double[i];

	// Berechne Konditionszahl und Fehler
	for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
	    A[i][j] = A_double[i][j];
	DenseMatrix<number> Ainv(A);
	inverse(Ainv);
	DenseMatrix<number> I(n,n);
	I.mm(Ainv,A);
	for (int i=0; i<n; i++) I[i][i]-= number(1.0);
	std::cout << "n=" << n
		  << " gamma=" << gamma
		  << " cond=" << A.norm_infty()*Ainv.norm_infty()
		  << " error=" << I.norm_infty()
		  << std::endl;
	conditionnumbers.push_back(A.norm_infty()*Ainv.norm_infty());
	errors.push_back(I.norm_infty());
      }
  std::cout << "condition numbers ";
  for (auto c : conditionnumbers) std::cout << " & " << c;
  std::cout << std::endl;
  std::cout << "errors ";
  for (auto e : errors) std::cout << " & " << e;
  std::cout << std::endl;
  
}
