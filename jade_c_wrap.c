// reference C-interfaces routines from 'jade_c.cxx' (libjade.a)
extern void wrap_jade (
		  double *B,	// write-only: B[nbc][nbc]
	      double *X,    // read-write: X[nbc][nbs]
	const int   nbc,
	const int   nbs   );

// pure C-interfaces calling into libjade.a
void jade (
		  double *icajade,	// write-only: icajade[dim][dim]
	      double *icacoffs,    // read-write: icacoffs[dim][num_samp]
	const int   dim,
	const int   num_samp  )
{
	wrap_jade(icajade, icacoffs, dim, num_samp);
}
