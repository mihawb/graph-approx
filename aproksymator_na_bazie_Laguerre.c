#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <float.h>

/*Oblicza wartosc odpowieniej funkcji Laguerra dla x*/
double fi(int n, double x)
{
    if(n == 0) return 1;
    else if(n == 1) return 1-x;
    else
    {
        int k = n - 1;
        return (((2*k+1-x)*fi(k,x))-k*fi(k-1,x))/(k+1);
    }
}

/*Pierwsza Pochodna*/
double dfi(int n, double x)
{
    double delta = 1.0e-6;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = fi(n,x1);
    double y2 = fi(n,x2);
    return (y2 - y1) / (x2 - x1);
}
/*Druga Pochodna*/
double d2fi(int n, double x)
{
    double delta = 1.0e-6;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = dfi(n,x1);
    double y2 = dfi(n,x2);
    return (y2 - y1) / (x2 - x1);
}

/*Trzecia Pochodna - v1.1*/
double d3fi(int n, double x)
{
    double delta = 1.0e-6;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = d2fi(n,x1);
    double y2 = d2fi(n,x2);
    return (y2 - y1) / (x2 - x1);
}


void make_spl(points_t * pts, spline_t * spl)
{
	matrix_t       *eqs= NULL;
        double         *x = pts->x;
        double         *y = pts->y;
        double          a = x[0];
        double          b = x[pts->n - 1];
        int             i, j, k;
        int             nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  	char *nbEnv= getenv( "APPROX_BASE_SIZE" );

        if( nbEnv != NULL && atoi( nbEnv ) > 0 )
                nb = atoi( nbEnv );

        eqs = make_matrix(nb, nb + 1);


	for (j = 0; j < nb; j++)
       	{
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(i, x[k]) * fi(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(j, x[k]));
	}
	
	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) 
			{
				double	ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (k, xx);
				spl->f1[i] += ck * dfi (k, xx);
				spl->f2[i] += ck * d2fi(k, xx);
				spl->f3[i] += ck * d3fi(k, xx);

			}
		}
	}
}
