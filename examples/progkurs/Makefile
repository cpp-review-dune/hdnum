include ../../make.def
HDNUMPATH  = ../..

# rule to build all programs without GMP support. That is the default
nogmp: ausgabe classes eingabe funktion gaussseidel hallohdnum matrizen mehrzahlen operationen pendel pendelmitfunktion pendelmitfunktionstemplate pendelmittimer pendelnumerisch pendelwhile string vektoren whitespace zahlen

all: nogmp

ausgabe: ausgabe.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

classes: classes.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

eingabe: eingabe.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

funktion: funktion.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

gaussseidel: gaussseidel.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

hallohdnum: hallohdnum.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

matrizen: matrizen.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

mehrzahlen: mehrzahlen.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

operationen: operationen.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

pendel: pendel.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

pendelmitfunktion: pendelmitfunktion.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

pendelmitfunktionstemplate: pendelmitfunktionstemplate.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

pendelmittimer: pendelmittimer.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

pendelnumerisch: pendelnumerisch.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

pendelwhile: pendelwhile.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

string: string.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

whitespace: whitespace.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

vektoren: vektoren.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

zahlen: zahlen.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

# clean up directory
clean:
	rm -f *.o ausgabe classes eingabe funktion gaussseidel hallohdnum matrizen mehrzahlen operationen pendel pendelmitfunktion pendelmitfunktionstemplate pendelmittimer pendelnumerisch pendelwhile string vektoren whitespace zahlen
