CC          = g++-mp-8
HDNUMPATH   = ../..
CCOPTS      = -I$(HDNUMPATH) -O3
#CCOPTS      = -I/opt/local/include -I$(HDNUMPATH) -O3
LIOPTSGMP   = -lgmpxx -lgmp -lm
#LIOPTSGMP   = -L/opt/local/lib -lgmpxx -lgmp -lm
LIOPTS      = -lm

corona_MAIN = corona
corona_SOURCES =

precision_MAIN = precision
precision_SOURCES =

fibonacci_MAIN = fibonacci
fibonacci_SOURCES =

vektoren_MAIN = vektoren
vektoren_SOURCES =

matrizen_MAIN = matrizen
matrizen_SOURCES =

lr_opcount_MAIN = lr_opcount
lr_opcount_SOURCES =

lr_MAIN = lr
lr_SOURCES =

wurzel_MAIN = wurzel
wurzel_SOURCES =

wurzelbanach_MAIN = wurzelbanach
wurzelbanach_SOURCES =

all : nogmp gmp

nogmp : $(corona_MAIN) $(precision_MAIN) $(vektoren_MAIN) $(matrizen_MAIN) $(fibonacci_MAIN) $(lr_opcount_MAIN)

gmp : $(lr_MAIN) $(wurzel_MAIN) $(wurzelbanach_MAIN)

$(zahlen_MAIN): $(zahlen_MAIN).o $(zahlen_SOURCES) Makefile
	$(CC) -o $(zahlen_MAIN) $(zahlen_MAIN).o $(zahlen_SOURCES) $(LIOPTS)

$(string_MAIN): $(string_MAIN).o $(string_SOURCES) Makefile
	$(CC) -o $(string_MAIN) $(string_MAIN).o $(string_SOURCES) $(LIOPTS)

$(pendelwhile_MAIN): $(pendelwhile_MAIN).o $(pendelwhile_SOURCES) Makefile
	$(CC) -o $(pendelwhile_MAIN) $(pendelwhile_MAIN).o $(pendelwhile_SOURCES) $(LIOPTS)

$(pendel_MAIN): $(pendel_MAIN).o $(pendel_SOURCES) Makefile
	$(CC) -o $(pendel_MAIN) $(pendel_MAIN).o $(pendel_SOURCES) $(LIOPTS)

$(pendelmittimer_MAIN): $(pendelmittimer_MAIN).o $(pendelmittimer_SOURCES) Makefile
	$(CC) -o $(pendelmittimer_MAIN) $(pendelmittimer_MAIN).o $(pendelmittimer_SOURCES) $(LIOPTS)

$(pendelnumerisch_MAIN): $(pendelnumerisch_MAIN).o $(pendelnumerisch_SOURCES) Makefile
	$(CC) -o $(pendelnumerisch_MAIN) $(pendelnumerisch_MAIN).o $(pendelnumerisch_SOURCES) $(LIOPTS)

$(pendelmitfunktion_MAIN): $(pendelmitfunktion_MAIN).o $(pendelmitfunktion_SOURCES) Makefile
	$(CC) -o $(pendelmitfunktion_MAIN) $(pendelmitfunktion_MAIN).o $(pendelmitfunktion_SOURCES) $(LIOPTS)

$(pendelmitfunktionstemplate_MAIN): $(pendelmitfunktionstemplate_MAIN).o $(pendelmitfunktionstemplate_SOURCES) Makefile
	$(CC) -o $(pendelmitfunktionstemplate_MAIN) $(pendelmitfunktionstemplate_MAIN).o $(pendelmitfunktionstemplate_SOURCES) $(LIOPTS)

$(whitespace_MAIN): $(whitespace_MAIN).o $(whitespace_SOURCES) Makefile
	$(CC) -o $(whitespace_MAIN) $(whitespace_MAIN).o $(whitespace_SOURCES) $(LIOPTS)

$(funktion_MAIN): $(funktion_MAIN).o $(funktion_SOURCES) Makefile
	$(CC) -o $(funktion_MAIN) $(funktion_MAIN).o $(funktion_SOURCES) $(LIOPTS)

$(corona_MAIN): $(corona_MAIN).o $(corona_SOURCES) Makefile
	$(CC) -o $(corona_MAIN) $(corona_MAIN).o $(corona_SOURCES) $(LIOPTS)

$(precision_MAIN): $(precision_MAIN).o $(precision_SOURCES) Makefile
	$(CC) -o $(precision_MAIN) $(precision_MAIN).o $(precision_SOURCES) $(LIOPTS)

$(hallohdnum_MAIN): $(hallohdnum_MAIN).o $(hallohdnum_SOURCES) Makefile
	$(CC) -o $(hallohdnum_MAIN) $(hallohdnum_MAIN).o $(hallohdnum_SOURCES) $(LIOPTS)

$(vektoren_MAIN): $(vektoren_MAIN).o $(vektoren_SOURCES) Makefile
	$(CC) -o $(vektoren_MAIN) $(vektoren_MAIN).o $(vektoren_SOURCES) $(LIOPTS)

$(matrizen_MAIN): $(matrizen_MAIN).o $(matrizen_SOURCES) Makefile
	$(CC) -o $(matrizen_MAIN) $(matrizen_MAIN).o $(matrizen_SOURCES) $(LIOPTS)

$(lr_opcount_MAIN): $(lr_opcount_MAIN).o $(lr_opcount_SOURCES) Makefile
	$(CC) -o $(lr_opcount_MAIN) $(lr_opcount_MAIN).o $(lr_opcount_SOURCES) $(LIOPTS)

$(modelproblem_MAIN): $(modelproblem_MAIN).o $(modelproblem_SOURCES) Makefile
	$(CC) -o $(modelproblem_MAIN) $(modelproblem_MAIN).o $(modelproblem_SOURCES) $(LIOPTS)

$(fibonacci_MAIN): $(fibonacci_MAIN).o $(fibonacci_SOURCES) Makefile
	$(CC) -o $(fibonacci_MAIN) $(fibonacci_MAIN).o $(fibonacci_SOURCES) $(LIOPTS)

$(ordertest_MAIN): $(ordertest_MAIN).o $(ordertest_SOURCES) Makefile
	$(CC) -o $(ordertest_MAIN) $(ordertest_MAIN).o $(ordertest_SOURCES) $(LIOPTS)


# Applications with GNU multiprecision library

$(eingabe_MAIN): $(eingabe_MAIN).o $(eingabe_SOURCES) Makefile
	$(CC) -o $(eingabe_MAIN) $(eingabe_MAIN).o $(eingabe_SOURCES) $(LIOPTSGMP)

$(mehrzahlen_MAIN): $(mehrzahlen_MAIN).o $(mehrzahlen_SOURCES) Makefile
	$(CC) -o $(mehrzahlen_MAIN) $(mehrzahlen_MAIN).o $(mehrzahlen_SOURCES) $(LIOPTSGMP)

$(lr_MAIN): $(lr_MAIN).o $(lr_SOURCES) Makefile
	$(CC) -o $(lr_MAIN) $(lr_MAIN).o $(lr_SOURCES) $(LIOPTSGMP)

$(wurzel_MAIN): $(wurzel_MAIN).o $(wurzel_SOURCES) Makefile
	$(CC) -o $(wurzel_MAIN) $(wurzel_MAIN).o $(wurzel_SOURCES) $(LIOPTSGMP)

$(wurzelbanach_MAIN): $(wurzelbanach_MAIN).o $(wurzelbanach_SOURCES) Makefile
	$(CC) -o $(wurzelbanach_MAIN) $(wurzelbanach_MAIN).o $(wurzelbanch_SOURCES) $(LIOPTSGMP)

.cc.o:
	$(CC) $(CCOPTS) -c $<

clean:
	rm -f *.o *.dat $(fibonacci_MAIN) $(precision_MAIN) $(vektoren_MAIN) $(matrizen_MAIN) $(lr_MAIN) $(wurzel_MAIN) $(wurzelbanach_MAIN)
