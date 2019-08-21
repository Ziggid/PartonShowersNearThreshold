CC=g++
CFLAGS=-c -Wall

fsshower : fsshower.o fsparton.o rng.o vec4d.o bindat.o 
	$(CC) fshower.o fsparton.o rng.o vec4d.o bindat.o -o fsshower

fsshower.o : shower.cpp
	$(CC) $(CFLAGS) shower.cpp -o fsshower.o

fsparton.o : partons.cpp
	$(CC) $(CFLAGS) partons.cpp -o fsparton.o

rng.o : rng.cpp
	$(CC) $(CFLAGS) rng.cpp

vec4d.o : vec4d.cpp
	$(CC) $(CFLAGS) vec4d.cpp

bindat.o : bindat.cpp
	$(CC) $(CFLAGS) bindat.cpp

clean:
	rm -rf *.o fsshower
