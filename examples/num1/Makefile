include ../../make.def
HDNUMPATH  = ../..

# rule to build all programs without GMP support. That is the default
nogmp: ecke eemodelproblem fem1d_linear fem1d_quadratic hodgkinhuxley iemodelproblem linearoscillator\
       modelproblem modelproblem_ordertest modelproblem_runge_kutta nbody ordertest system_example vanderpol

# rule to build programs with GMP support
gmp: modelproblem_high_dim

all: nogmp gmp

#
# Programs reqiring GMP support
#
ecke: ecke.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

eemodelproblem: eemodelproblem.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

fem1d_linear: fem1d_linear.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

fem1d_quadratic: fem1d_quadratic.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

hodgkinhuxley: hodgkinhuxley.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

iemodelproblem: iemodelproblem.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

linearoscillator: linearoscillator.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

modelproblem: modelproblem.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

modelproblem_ordertest: modelproblem_ordertest.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

modelproblem_runge_kutta: modelproblem_runge_kutta.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

nbody: nbody.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

ordertest: ordertest.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

system_example: system_example.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

vanderpol: vanderpol.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

#
# Programs reqiring GMP support
#

modelproblem_high_dim: modelproblem_high_dim.cc
	$(CC) $(CCFLAGS) $(GMPCCFLAGS) -o $@ $^ $(GMPLFLAGS) $(LFLAGS)

# clean up directory
clean:
	rm -f *.o ecke eemodelproblem fem1d_linear fem1d_quadratic hodgkinhuxley iemodelproblem linearoscillator\
       modelproblem modelproblem_ordertest modelproblem_runge_kutta nbody ordertest system_example vanderpol modelproblem_high_dim

