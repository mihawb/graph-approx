#include "makespl.h"
#include "gaus/piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
algorytm dla Wielomianów Laguerre’a
Lk(x) = ((2k-1-x)/k) * Lk-1(x) - ((k-1)/k) * Lk-2(x).

cz1 => ((2k-1-x)/k) * Lk-1(x)
cz2 => ((k - 1)/k)Lk-2(x)


jako że bez aproxymacji niepotrzebne nam  (n - liczba funkcji , a,b - granice przedzialu aproksymacji) ale zostawione aby zgodny interfejs został

   a,b - granice przedzialu aproksymacji
   n - liczba funkcji
   i  - numer funkcji
   x - wartość dla ktorej obliczana jest wartosc funkcji
*/
double fi(double a, double b, int n, int i, double x) {
	if (i == 0)
		return 1;

	if (i == 1)
		return (-x + 1);


	double cz1 = (double)(2 * i - 1 - x) / i;
	cz1 = cz1 * fi(a, b, n, i - 1, x);  /* wywolanie rekurencyjne dla k-1	*/

	double cz2 = (double)(i - 1) / i;
	cz2 = cz2 * fi(a, b, n, i - 2, x); /* wywolanie rekurencyjne dla k-2	*/

	double result = cz1 - cz2;

	return result;

}

/*
numeryczna aproxymacja pierwszej pochodnej przy pomocy centralnej 5 punktowej:
 progresywne:
	!2 punktowej  ==>  f’(x) = [f(x+h) - f(x) ]/h
	 3 punktowej  ==>  f’(x) = [f(x+h) - f(x-h)]/2h
	 5-punktowej ==> f’(x) = [ f(x-2h) - 8f(x-h) +8f(x+h) -f(x+2h)]/12h

centralne:
	f’(x) = (-f(x+2h) + 8*f(x+h) - 8*f(x-1) + f(x-2h)) / (12 * h)


   a,b - granice przedzialu aproksymacji
   n - liczba funkcji
   i  - numer funkcji
   x - wartość dla ktorej obliczana jest wartosc funkcji
*/

double dfi(double a, double b, int n, int i, double x) {     //5 punktowa centralna
	/*przedzial*/
	double h = (b - a) / (n - 1);

	/*nastepny x*/
	double x_next_1 = x + h;
	double x_next_2 = x + h + h;

	/*poprzedni x*/
	double x_prev_1 = x - h;
	double x_prev_2 = x - h - h;


	return (-fi(a, b, n, i, x_next_2) + 8 * fi(a, b, n, i, x_next_1) - 8 * fi(a, b, n, i, x_prev_1) + fi(a, b, n, i, x_prev_2)) / (12 * h);
}

/* 2 punktowa progresywna obardczona duzym bledem slabe wykresy (nie jest uzyta)*/
double dfi_2p(double a, double b, int n, int i, double x) {    //2 punktowa progresywna
	/*przedzial*/
	double h = (b - a) / (n - 1);

	/*nastepny x*/
	double x1 = x + h;

	return (fi(a, b, n, i, x1) - fi(a, b, n, i, x)) / h;
}


/*
 numeryczna aproxymacja drógiej pochodnej przy pomocy aproksymacji 5 punktowej centralnej:
	f''(x) = -f(x+2h) + 16*f(x+1h) - 30*f(x) + 16*f(x-1h) - f(x-2h)) / (12*h*h)

   n - liczba funkcji
   a,b - granice przedzialu aproksymacji
   i  - numer funkcji
   x - wspolrzedna dla ktorej obliczana jest wartosc funkcji

*/

double d2fi(double a, double b, int n, int i, double x) {     // 5 punktowa
	/* można tych zmiennych nie delkarować ale dla czytelności zostły*/

	/*przedzial*/
	double h = (b - a) / (n - 1);

	/*nastepny x*/
	double x_next_1 = x + h;
	double x_next_2 = x + h + h;

	/*poprzedni x*/
	double x_prev_1 = x - h;
	double x_prev_2 = x - h - h;

	return (-fi(a, b, n, i, x_next_2) + 16 * fi(a, b, n, i, x_next_1) - 30 * fi(a, b, n, i, x) + 16 * fi(a, b, n, i, x_prev_1) - fi(a, b, n, i, x_prev_2)) / (12 * h * h);
}

/* 3 punktowa progresywna obardczona duzym bledem slabe wykresy (nie jest uzyta)*/
double d2fi_3p(double a, double b, int n, int i, double x) {
	/* można tych zmiennych nie delkarować ale dla czytelności zostły*/

	/*przedzial*/
	double h = (b - a) / (n - 1);

	/*nastepny x*/
	double x1 = x + h;

	/*poprzedni x*/
	double x2 = x + h +h;

	return (fi(a, b, n, i, x2) - (2 * fi(a, b, n, i, x1)) + fi(a, b, n, i, x)) / (h * h);
}

/*
 numeryczna aproxymacja trzeciej pochodnej przy pomocy aproksymacji 4 punktowej centralnej+progres:
	- f(x+3h) + 8*f(x+2h) - 13*f(x+1h) + 13*f(x-1h) - 8*f(x-2h) + f(x-3h)) / (8*h*h*h)

   n - liczba funkcji
   a,b - granice przedzialu aproksymacji
   i  - numer funkcji
   x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
*/

double d3fi(double a, double b, int n, int i, double x) {     //6 punktowa
	/*przedzial*/
	double h = (b - a) / (n - 1);

	/*nastepny x*/
	double x1 = x + h;

	/*nastepny x*/
	double x_next_1 = x + h;
	double x_next_2 = x + h + h;
	double x_next_3 = x + h + h + h;

	/*poprzedni x*/
	double x_prev_1 = x - h;
	double x_prev_2 = x - h - h;
	double x_prev_3 = x - h - h - h;

	return ( - fi(a, b, n, i, x_next_3) + 8 * fi(a, b, n, i, x_next_2) - 13 * fi(a, b, n, i, x_next_1) + 13 * fi(a, b, n, i, x_prev_1) - 8 * fi(a, b, n, i, x_prev_2) + fi(a, b, n, i, x_prev_3)) / (8 * h * h * h);
}
 
/* 4 punktowa progresywna obardczona duzym bledem slabe wykresy (nie jest uzyta)*/
double d3fi_4p(double a, double b, int n, int i, double x) {
	/*przedzial*/
	double h = (b - a) / (n - 1);

	/*nastepny x*/
	double x1 = x + h;

	/*nastepny x +h+h*/
	double x2 = x + h + h;

	/*poprzedni x*/
	double x3 = x + h +h +h;


	return (fi(a, b, n, i, x3) - (3 * fi(a, b, n, i, x2)) + (3 * fi(a, b, n, i, x1)) - fi(a, b, n, i, x)) / (h * h * h);
}


/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %d", hi[j] );
	fprintf( out, "] hx=[" );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %g", hx[j] );
	fprintf( out, "]\n" );
}

void make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for( j= 0; j < nb; j++ )
			xfi( a, b, nb, j, tst );
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
				fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, nb) == 0) {
	
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			double yi= 0;
			double dyi= 0;
			double d2yi= 0;
			double d3yi= 0;
			double xi= a + i * dx;
			for( k= 0; k < nb; k++ ) {
							yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}
