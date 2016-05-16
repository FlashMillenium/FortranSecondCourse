CC=gfortran
CFLAGS=-Wall -std=f2008ts -Wno-maybe-uninitialized -flto

all:
	$(CC) $(CFLAGS) -c src/FEHL.F -J obj/ -o obj/feh.o
	$(CC) $(CFLAGS) -c src/RKFS.F -J obj/ -o obj/rek.o
	$(CC) $(CFLAGS) -c src/RKF45.F -J obj/ -o obj/dec.o
	$(CC) $(CFLAGS) -c src/main.f95 -I obj/ -o obj/test.o
	$(CC) $(CFLAGS) -o bin/difffound ./obj/feh.o obj/rek.o obj/dec.o obj/test.o
clean:
	rm -f bin/*.*
	rm -f obj/*.*
	rm -f src/*.swp
	rm -f src/*.swm
	rm -f src/*.swn
run:
	cd ./bin; ./difffound
	
