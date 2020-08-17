// operationen.cc
#include <iostream>

int main ()
{
  int a = 6; //Variable a vom Typ int wird definiert.
  int b = 2;
  int c = a + b; // Wert von a + b wird c zugewiesen.
	
  // Natürlich geht auch:
  c = a*b;  //c wurde oben bereits definiert
  c = a-b;
  c = a/b;  // nur möglich, falls a/b ganzzahlig sonst wird gerundet!

  // Auch möglich:
  c = a/c;  //a wird durch c geteilt und dann c zugewiesen.

  // Oder eben auch komplexere Ausdrücke:
  c = b+a/b-a;
  c = (a+b)/(a-b);
  // es gilt Punkt-vor-Strich, aber Klammern gehen vor

  return 0;
}
