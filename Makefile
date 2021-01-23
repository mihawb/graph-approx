aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

aproxl: main.o splines.o points.o aproksymator_na_bazie_laguerra2.o gaus/libge.a
	$(CC) -o aproxl main.o splines.o points.o aproksymator_na_bazie_laguerra2.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

aproksymator_na_bazie_laguerra2.o: makespl.h points.h gaus/piv_ge_solver.h aproksymator_na_bazie_laguerra2.c
	$(CC) -I gaus -c aproksymator_na_bazie_laguerra2.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

.PHONY: clean test

clean:
	-rm *.o aprox intrp prosta aproxl
	-rm do*plot spldo*

test:
	./aproxl -s spldo2 -p test/dane.do2 -g do2plot -f -3 -t 3 -n 500
	./aproxl -s spldo3 -p test/dane.do3 -g do3plot -f -3 -t 3 -n 500
	./aproxl -s spldo4 -p test/dane.do4 -g do4plot -f -3 -t 3 -n 500
	./aproxl -s spldo5 -p test/dane.do5 -g do5plot -f -3 -t 3 -n 500