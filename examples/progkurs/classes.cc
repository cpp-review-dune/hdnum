// classes.cc
#include <iostream> // header für Ein-/Ausgabe
#include <cmath>    // mathematische Funktionen

class Vector3d // einfaches Beispiel für eine Klasse
{
public:
  float x; //Klassenvariablen
  float y;
  float z;

  Vector3d(float par1, float par2, float par3)
  {
    x = par1;
    y = par2;
    z = par3;
  }

  float getNorm() //Klassenmethode (Norm des Vektors)
  {
    return sqrt(x*x + y*y + z*z);
  }

  float dot(Vector3d U) // Skalarprodukt mit anderem Vektor
  {
    return x*U.x + y*U.y + z*U.z;
  }
};


int main ()
{
  Vector3d V(1.0,2.0,3.0); // Erstellen einer Variable der Klasse Vector3d
  Vector3d W(-1.0,-1.0,1.0); // Erstellen einer anderen Variable der Klasse Vector3d
  std::cout << V.x << std::endl; // Zugriff auf Variable x von V
  std::cout << W.x << std::endl; // Zugriff auf Variable x von W
  std::cout << V.getNorm() << std::endl; // Norm von V
  std::cout << W.getNorm() << std::endl; // Norm von W
  std::cout << V.dot(W) << std::endl;    // Skalarprodukt zwischen V und W
  std::cout << W.dot(V) << std::endl;    // Skalarprodukt zwischen W und V (das selbe)

  return 0;
}
