CC          = g++
HDNUMPATH   = ../..
CCOPTS      = -I$(HDNUMPATH) -O3
#CCOPTS      = -I/opt/local/include -I$(HDNUMPATH) -O3
LIOPTSGMP   = -lgmpxx -lgmp -lm
#LIOPTSGMP   = -L/opt/local/lib -lgmpxx -lgmp -lm

modelproblem_MAIN = modelproblem
modelproblem_SOURCES =

modelproblem_high_dim_MAIN = modelproblem_high_dim
modelproblem_high_dim_SOURCES =

modelproblem_ordertest_MAIN = modelproblem_ordertest
modelproblem_ordertest_SOURCES =

modelproblem_runge_kutta_MAIN = modelproblem_runge_kutta
modelproblem_runge_kutta_SOURCES =

system_example_MAIN = system_example
system_example_SOURCES =

ordertest_MAIN = ordertest
ordertest_SOURCES =

hodgkinhuxley_MAIN = hodgkinhuxley
hodgkinhuxley_SOURCES =

nbody_MAIN = nbody
nbody_SOURCES =

vanderpol_MAIN = vanderpol
vanderpol_SOURCES =

wurzel_MAIN = wurzel
wurzel_SOURCES =

iemodelproblem_MAIN = iemodelproblem
iemodelproblem_SOURCES =

ecke_MAIN = ecke
ecke_SOURCES =

all : nogmp gmp

nogmp : $(modelproblem_MAIN) $(modelproblem_ordertest_MAIN) $(modelproblem_runge_kutta_MAIN) $(system_example_MAIN) $(ordertest_MAIN) $(hodgkinhuxley_MAIN) $(nbody_MAIN) $(vanderpol_MAIN) $(iemodelproblem_MAIN) $(ecke_MAIN)

gmp : $(modelproblem_high_dim_MAIN)

$(modelproblem_MAIN): $(modelproblem_MAIN).o $(modelproblem_SOURCES) Makefile
	$(CC) -o $(modelproblem_MAIN) $(modelproblem_MAIN).o $(modelproblem_SOURCES) $(LIOPTS)

$(modelproblem_high_dim_MAIN): $(modelproblem_high_dim_MAIN).o $(modelproblem_high_dim_SOURCES) Makefile
	$(CC) -o $(modelproblem_high_dim_MAIN) $(modelproblem_high_dim_MAIN).o $(modelproblem_high_dim_SOURCES) $(LIOPTSGMP)

$(modelproblem_ordertest_MAIN): $(modelproblem_ordertest_MAIN).o $(modelproblem_ordertest_SOURCES) Makefile
	$(CC) -o $(modelproblem_ordertest_MAIN) $(modelproblem_ordertest_MAIN).o $(modelproblem_ordertest_SOURCES) $(LIOPTS)

$(modelproblem_runge_kutta_MAIN): $(modelproblem_runge_kutta_MAIN).o $(modelproblem_runge_kutta_SOURCES) Makefile
	$(CC) -o $(modelproblem_runge_kutta_MAIN) $(modelproblem_runge_kutta_MAIN).o $(modelproblem_runge_kutta_SOURCES) $(LIOPTS)

$(system_example_MAIN): $(system_example_MAIN).o $(system_example_SOURCES) Makefile
	$(CC) -o $(system_example_MAIN) $(system_example_MAIN).o $(system_example_SOURCES) $(LIOPTS)

$(ordertest_MAIN): $(ordertest_MAIN).o $(ordertest_SOURCES) Makefile
	$(CC) -o $(ordertest_MAIN) $(ordertest_MAIN).o $(ordertest_SOURCES) $(LIOPTS)

$(hodgkinhuxley_MAIN): $(hodgkinhuxley_MAIN).o $(hodgkinhuxley_SOURCES) Makefile
	$(CC) -o $(hodgkinhuxley_MAIN) $(hodgkinhuxley_MAIN).o $(hodgkinhuxley_SOURCES) $(LIOPTS)

$(nbody_MAIN): $(nbody_MAIN).o $(nbody_SOURCES) Makefile
	$(CC) -o $(nbody_MAIN) $(nbody_MAIN).o $(nbody_SOURCES) $(LIOPTS)

$(vanderpol_MAIN): $(vanderpol_MAIN).o $(vanderpol_SOURCES) Makefile
	$(CC) -o $(vanderpol_MAIN) $(vanderpol_MAIN).o $(vanderpol_SOURCES) $(LIOPTS)

$(iemodelproblem_MAIN): $(iemodelproblem_MAIN).o $(iemodelproblem_SOURCES) Makefile
	$(CC) -o $(iemodelproblem_MAIN) $(iemodelproblem_MAIN).o $(iemodelproblem_SOURCES) $(LIOPTS)

$(ecke_MAIN): $(ecke_MAIN).o $(ecke_SOURCES) Makefile
	$(CC) -o $(ecke_MAIN) $(ecke_MAIN).o $(ecke_SOURCES) $(LIOPTS)

# Applications with GNU multiprecision library

$(eingabe_MAIN): $(eingabe_MAIN).o $(eingabe_SOURCES) Makefile
	$(CC) -o $(eingabe_MAIN) $(eingabe_MAIN).o $(eingabe_SOURCES) $(LIOPTSGMP)

.cc.o:
	$(CC) $(CCOPTS) -c $<

clean:
	rm -f *.o *.dat *.gp $(modelproblem_MAIN) $(modelproblem_high_dim_MAIN) $(modelproblem_ordertest_MAIN) $(modelproblem_runge_kutta_MAIN) $(system_example_MAIN) $(ordertest_MAIN) $(hodgkinhuxley_MAIN) $(nbody_MAIN) $(vanderpol_MAIN) $(iemodelproblem_MAIN) $(ecke_MAIN)
