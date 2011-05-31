all: virus_gencor

virus_gencor: virus_gencor.o
	mpic++ -o virus_gencor virus_gencor.o

virus_gencor.o: virus_gencor.cc virus_gencor.h
	mpic++ -c virus_gencor.cc

clean:
	rm -f virus_gencor \
		virus_gencor.o
