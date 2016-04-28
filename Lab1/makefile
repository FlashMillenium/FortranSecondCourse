CC = gfortran
CFLAGS=-Wall -std=f2008ts -fimplicit-none -flto

all:
	$(CC) $(CFLAGS) -c src/SPLINE.f90 -J obj/ -o obj/SPLINE.o
	$(CC) $(CFLAGS) -c src/QUANC8.F -J obj/ -o obj/quan.o
	$(CC) $(CFLAGS) -c src/SEVAL.F -J obj/ -o obj/seval.o
	$(CC) $(CFLAGS) -c src/main.f95 -I obj/ -o obj/test.o
	$(CC) $(CFLAGS) -o bin/lab.darn ./obj/quan.o obj/seval.o obj/SPLINE.o obj/test.o 
clean:
	rm -f bin/*.*
	rm -f obj/*.*
run:
	cd ./bin; ./lab.darn
