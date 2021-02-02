#include <math.h>
#include "nrutil.h"
#include "nr.h"
#include <functional>
#define MAXITS 20000
#define EPS 1.0e-30
#define TOLF 1.0e-30
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-30
//#define float double
/* Here MAXITS is the maximum number of iterations; EPS is a number close to the machine
precision; TOLF is the convergence criterion on function values; TOLX is the convergence criterion
on δx; STPMX is the scaled maximum step length allowed in line searches; TOLMIN is used to
decide whether spurious convergence to a minimum of fmin has occurred. */
#define FREERETURN {free_dvector(fvec,1,n);free_dvector(xold,1,n);\
    free_dvector(w,1,n);free_dvector(t,1,n);free_dvector(s,1,n);\
    free_dmatrix(r,1,n,1,n);free_dmatrix(qt,1,n,1,n);free_dvector(p,1,n);\
    free_dvector(g,1,n);free_dvector(fvcold,1,n);free_dvector(d,1,n);\
    free_dvector(c,1,n);return;}
    
int nn;                                                                                      // Global variables to communicate with fmin.
double *fvec;
//void (*nrfuncv)(int n, double v[], double f[]);

std::function<void(int, double [], double [])> nrfuncv;

void broydn(double x[], int n, int *check,
    std::function<void(int, double [], double [])> vecfunc
	//void (*vecfunc)(int, double [], double [])
	)
/* Given an initial guess x[1..n] for a root in n dimensions, find the root by Broyden’s method
embedded in a globally convergent strategy. The vector of functions to be zeroed, called
fvec[1..n] in the routine below, is returned by the user-supplied routine vecfunc(n,x,fvec) .
The routine fdjac and the function fmin from newt are used. The output quantity check
is false (0) on a normal return and true (1) if the routine has converged to a local minimum
of the function fmin or if Broyden’s method can make no further progress. In this case try
restarting from a different initial guess. */
{
    //void fdjac(int n, double x[], double fvec[], double **df,
	//std::function<void(int, double [], double [])> vecfunc);
    double fmin1(double x[]);
    void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
        double *f, double stpmax, int *check, double (*func)(double []));
    void qrdcmp(double **a, int n, double *c, double *d, int *sing);
    void qrupdt(double **r, double **qt, int n, double u[], double v[]);
    void rsolv(double **a, int n, double d[], double b[]);
    int i,its,j,k,restrt,sing,skip;
    double den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
    double *g,*p,**qt,**r,*s,*t,*w,*xold;
    
    c=dvector(1,n);
    d=dvector(1,n);
    fvcold=dvector(1,n);
    g=dvector(1,n);
    p=dvector(1,n);
    qt=dmatrix(1,n,1,n);
    r=dmatrix(1,n,1,n);
    s=dvector(1,n);
    t=dvector(1,n);
    w=dvector(1,n);
    xold=dvector(1,n);
    fvec=dvector(1,n);                                                    // Define global variables.
    nn=n;
    nrfuncv=vecfunc;
    f=fmin1(x);                                                                // The vector fvec is also computed by this
    test=0.0;                                                                   // call.
    for (i=1;i<=n;i++)                                                  // Test for initial guess being a root. Use more
        if (fabs(fvec[i]) > test)test=fabs(fvec[i]);   // stringent test than sim-
    if (test < 0.01*TOLF) {                                           // ply TOLF.
        *check=0;
        FREERETURN
    }
    for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);               // Calculate stpmax for line searches.
    stpmax=STPMX*FMAX(sqrt(sum),(double)n);
    restrt=1;                                                                                     // Ensure initial Jacobian gets computed.
    for (its=1;its<=MAXITS;its++) {                                           // Start of iteration loop.
        if (restrt) {
            fdjac(n,x,fvec,r,vecfunc);                                                 // Initialize or reinitialize Jacobian in r.
            qrdcmp(r,n,c,d,&sing);                                                      // QR decomposition of Jacobian.
            if (sing) nrerror("singular Jacobian in broydn");
            for (i=1;i<=n;i++) {                                                          // Form Q T explicitly.
                for (j=1;j<=n;j++) qt[i][j]=0.0;
                qt[i][i]=1.0;
            }
            for (k=1;k<n;k++) {
                if (c[k]) {
                    for (j=1;j<=n;j++) {
                        sum=0.0;
                        for (i=k;i<=n;i++)
                            sum += r[i][k]*qt[i][j];
                        sum /= c[k];
                        for (i=k;i<=n;i++)
                            qt[i][j] -= sum*r[i][k];
                    }
                }
            }
            for (i=1;i<=n;i++) {                                 // Form R explicitly.
                r[i][i]=d[i];
                for (j=1;j<i;j++) r[i][j]=0.0;
            }
        } else {                                                             // Carry out Broyden update.
            for (i=1;i<=n;i++) s[i]=x[i]-xold[i];         // s = δx.
            for (i=1;i<=n;i++) {                                      // t = R · s.
                for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
                t[i]=sum;
            }
            skip=1;
            for (i=1;i<=n;i++) {                                                            // w = δF − B · s.
                for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
                w[i]=fvec[i]-fvcold[i]-sum;
                if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) skip=0;  
                // Don’t update with noisy components of w.
                else w[i]=0.0;
            }
            if (!skip) {
                for (i=1;i<=n;i++) {                                                     // t = Q T · w.
                    for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
                    t[i]=sum;
                }
                for (den=0.0,i=1;i<=n;i++) den += SQR(s[i]);
                for (i=1;i<=n;i++) s[i] /= den;                                     // Store s/(s · s) in s.
                qrupdt(r,qt,n,t,s);                                                             // Update R and Q T .
                for (i=1;i<=n;i++) {
                    if (r[i][i] == 0.0) nrerror("r singular in broydn");
                    d[i]=r[i][i];                                                                    // Diagonal of R stored in d.
                }
            }
        }
        for (i=1;i<=n;i++) {                                                          // Compute ∇f ≈ (Q · R) T · F for the line search.
            for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
            g[i]=sum;
        }
        for (i=n;i>=1;i--) {
            for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
            g[i]=sum;
        }
        for (i=1;i<=n;i++) {                                                         // Store x and F.
            xold[i]=x[i];
            fvcold[i]=fvec[i];
        }
        fold=f;                                                                                  // Store f .
        for (i=1;i<=n;i++) {                                                         // Right-hand side for linear equations is −Q T · F.
            for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
            p[i] = -sum;
        }
        rsolv(r,n,d,p);                                                                // Solve linear equations.
        lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin1);
        // lnsrch returns new x and f . It also calculates fvec at the new x when it calls fmin.
        test=0.0;                                  // Test for convergence on function values.
        for (i=1;i<=n;i++)
            if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
        if (test < TOLF) {
            *check=0;
            FREERETURN
        }
        if (*check) {                                                // True if line search failed to find a new x.
            if (restrt) FREERETURN                        // Failure; already tried reinitializing the Jaco-
            else {                                                    // bian.
                test=0.0;                                          // Check for gradient of f zero, i.e., spurious
                den=FMAX(f,0.5*n);                      // convergence.
                for (i=1;i<=n;i++) {
                    temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
                    if (temp > test) test=temp;
                }
                if (test < TOLMIN) FREERETURN
                else restrt=1;                                    // Try reinitializing the Jacobian.
            }
        } else {                                      // Successful step; will use Broyden update for
            restrt=0;                                // next step.
            test=0.0;                                // Test for convergence on δx.
            for (i=1;i<=n;i++) {
                temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
                if (temp > test) test=temp;
            }
            if (test < TOLX) FREERETURN
        }
    }
    nrerror("MAXITS exceeded in broydn");
    FREERETURN
}

void rsolv(double **a, int n, double d[], double b[])
/* Solves the set of n linear equations R · x = b, where R is an upper triangular matrix stored in
     a and d . a[1..n][1..n] and d[1..n] are input as the output of the routine qrdcmp and
     are not modified. b[1..n] is input as the right-hand side vector, and is overwritten with the
     solution vector on output. */
{
    int i,j;
    double sum;
    
    b[n] /= d[n];
    for (i=n-1;i>=1;i--) {
        for (sum=0.0,j=i+1;j<=n;j++) sum += a[i][j]*b[j];
        b[i]=(b[i]-sum)/d[i];
    }
}

void qrdcmp(double **a, int n, double *c, double *d, int *sing)
/* Constructs the QR decomposition of a[1..n][1..n] . The upper triangular matrix R is re-
     turned in the upper triangle of a , except for the diagonal elements of R which are returned in
     d[1..n] . The orthogonal matrix Q is represented as a product of n − 1 Householder matrices
     Q 1 . . . Q n−1 , where Q j = 1 − u j ⊗ u j /c j . The ith component of u j is zero for i = 1, . . . , j − 1
     while the nonzero components are returned in a[i][j] for i = j, . . . , n. sing returns as
     true ( 1 ) if singularity is encountered during the decomposition, but the decomposition is still
     completed in this case; otherwise it returns false ( 0 ). */
{
    int i,j,k;
    double scale,sigma,sum,tau;

    *sing=0;
    for (k=1;k<n;k++) {
        scale=0.0;
        for (i=k;i<=n;i++) scale=FMAX(scale,fabs(a[i][k]));
        if (scale == 0.0) {                    // Singular case.
            *sing=1;
            c[k]=d[k]=0.0;
        } else {                                       // Form Q k and Q k · A.
            for (i=k;i<=n;i++) a[i][k] /= scale;
            for (sum=0.0,i=k;i<=n;i++) sum += SQR(a[i][k]);
            sigma=SIGN(sqrt(sum),a[k][k]);
            a[k][k] += sigma;
            c[k]=sigma*a[k][k];
            d[k] = -scale*sigma;
            for (j=k+1;j<=n;j++) {
                for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j];
                tau=sum/c[k];
                for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
            }
        }
    }
    d[n]=a[n][n];
    if (d[n] == 0.0) *sing=1;
}

void qrupdt(double **r, double **qt, int n, double u[], double v[])
/* Given the QR decomposition of some n × n matrix, calculates the QR decomposition of the
     matrix Q · (R+ u ⊗ v). The quantities are dimensioned as r[1..n][1..n] , qt[1..n][1..n] ,
     u[1..n] , and v[1..n] . Note that Q T is input and returned in qt. */
{
    void rotate(double **r, double **qt, int n, int i, double a, double b);
    int i,j,k;
    for (k=n;k>=1;k--) {                            // Find largest k such that u[k]  = 0.
        if (u[k]) break;
    }
    if (k < 1) k=1;
    for (i=k-1;i>=1;i--) {                            // Transform R + u ⊗ v to upper Hessenberg.
        rotate(r,qt,n,i,u[i],-u[i+1]);
        if (u[i] == 0.0) u[i]=fabs(u[i+1]);
        else if (fabs(u[i]) > fabs(u[i+1]))
            u[i]=fabs(u[i])*sqrt(1.0+SQR(u[i+1]/u[i]));
        else u[i]=fabs(u[i+1])*sqrt(1.0+SQR(u[i]/u[i+1]));
    }
    for (j=1;j<=n;j++) r[1][j] += u[1]*v[j];
    for (i=1;i<k;i++)                                    // Transform upper Hessenberg matrix to upper tri-
        rotate(r,qt,n,i,r[i][i],-r[i+1][i]);    // angular.
}

#include <math.h>
#include "nrutil.h"
void rotate(double **r, double **qt, int n, int i, double a, double b)
/* Given matrices r[1..n][1..n] and qt[1..n][1..n] , carry out a Jacobi rotation
     √ on rows 
     i and i + 1 √ of each matrix. a and b are the parameters of the rotation: cos θ = a/ a 2 + b 2 ,
     sin θ = b/ a 2 + b 2. */
{
    int j;
    double c,fact,s,w,y;
    if (a == 0.0) {                             // Avoid unnecessary overflow or underflow.
        c=0.0;
        s=(b >= 0.0 ? 1.0 : -1.0);
    } else if (fabs(a) > fabs(b)) {
        fact=b/a;
        c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
        s=fact*c;
    } else {
        fact=a/b;
        s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
        c=fact*s;
    }
    for (j=i;j<=n;j++) {                            // Premultiply r by Jacobi rotation.
        y=r[i][j];
        w=r[i+1][j];
        r[i][j]=c*y-s*w;
        r[i+1][j]=s*y+c*w;
    }
    for (j=1;j<=n;j++) {                           // Premultiply qt by Jacobi rotation.
        y=qt[i][j];
        w=qt[i+1][j];
        qt[i][j]=c*y-s*w;
        qt[i+1][j]=s*y+c*w;
    }
}