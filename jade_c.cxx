// C++ Templates of Problem Set #1 subroutines
//
// Generates 'libjade_c.a'
#include "jade_c.h"
#include <cstdlib>
#include <iostream>

// C-interface routines for F2PY wrapper

extern "C" void wrap_jade (
		  double *B,	// write-only: B[nbc][nbc]
	      double *X,    // read-write: X[nbc][nbs]
	const int   nbc,
	const int   nbs)
{
	try {
		(void)ce::jade(B, X, nbc, nbs);
	} catch(const char *e) {
		std::cerr << "jade: " << e << std::endl;
		std::exit(1);
	} catch(...) {
		std::cerr << "jade: Unhandled Exception" << std::endl;
		std::exit(1);
	}
}
