include ../make.def
HDNUMPATH  = ..

# rule to build all programs without GMP support. That is the default
nogmp: example

# rule to build programs with GMP support
gmp: wurzel

all: nogmp gmp

#
# Programs reqiring no GMP support
#

example: example.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(LFLAGS)

#
# Programs reqiring GMP support
# just an example to see how it would work, there is no wurzel.cc here
#

wurzel: wurzel.cc
	$(CC) $(CCFLAGS) $(GMPCCFLAGS) -o $@ $^ $(GMPLFLAGS) $(LFLAGS)

# clean up directory
clean:
	rm -f *.o example
