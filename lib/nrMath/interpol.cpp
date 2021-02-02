#include "nr.h"
#include <math.h>

#define TINY 1.0e-25                // A small number.
#define FREERETURN {free_vector(d,1,n);free_vector(c,1,n);return;}
#define float double
#include "nrutil.h"

void ratint(float xa[], float ya[], int n, float x, float *y, float *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of
y and an accuracy estimate dy. The value returned is that of the diagonal rational function,
evaluated at x, which passes through the n points (xa_i,ya_i), i=1...n. */
{
	int m,i,ns=1;
	float w,t,hh,h,dd,*c,*d;
	
	c=dvector(1,n);
	d=dvector(1,n);
	hh=fabs(x-xa[1]);
	for (i=1;i<=n;i++) {
		h=fabs(x-xa[i]);
		if (h == 0.0) {
			*y=ya[i];
			*dy=0.0;
			FREERETURN
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;            // The TINY part is needed to prevent a rare zero-over-zero
	}                               // condition.
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;            // h will never be zero, since this was tested in the initial-
			t=(xa[i]-x)*d[i]/h;     // izing loop.
			dd=t-c[i+1];
			if (dd == 0.0) nrerror("Error in routine ratint");
			/* This error condition indicates that the interpolating function has a pole at the
			requested value of x. */
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	FREERETURN
}