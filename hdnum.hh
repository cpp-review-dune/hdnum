// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_HDNUM_HH
#define HDNUM_HDNUM_HH

#define HDNUM_DEVEL_MODE 1

#include <math.h>

#include <complex>

#if HDNUM_HAS_GMP
#include <gmpxx.h>
#endif

#include "src/densematrix.hh"
#include "src/exceptions.hh"
#include "src/lr.hh"
#include "src/newton.hh"
#include "src/ode.hh"
#include "src/opcounter.hh"
#include "src/pde.hh"
#include "src/precision.hh"
#include "src/qr.hh"
#include "src/rungekutta.hh"
#include "src/sgrid.hh"
#include "src/sparsematrix.hh"
#include "src/timer.hh"
#include "src/vector.hh"

#endif
