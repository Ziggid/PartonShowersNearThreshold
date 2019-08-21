CC=g++
CFLAGS=-c -Wall

isshower : isshower.o isparton.o rng.o ConfigFile.o vec4.o bindat.o mstwpdf.o
	$(CC) isshower.o isparton.o rng.o ConfigFile.o vec4.o bindat.o mstwpdf.o -o isshower

isshower.o : isshower.cpp
	$(CC) $(CFLAGS) isshower.cpp

isparton.o : isparton.cpp
	$(CC) $(CFLAGS) isparton.cpp

rng.o : rng.cpp
	$(CC) $(CFLAGS) rng.cpp

ConfigFile.o : ConfigFile.cpp
	$(CC) $(CFLAGS) ConfigFile.cpp

vec4.o : vec4.cpp
	$(CC) $(CFLAGS) vec4.cpp

bindat.o : bindat.cpp
	$(CC) $(CFLAGS) bindat.cpp

mstwpdf.o : mstwpdf.cc
	$(CC) $(CFLAGS) mstwpdf.cc

clean:
	rm -rf *.o isshower
