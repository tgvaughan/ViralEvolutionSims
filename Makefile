all: virus_gencor_simple

tags: virus_gencor_simple.cc poissonian.cc poissonian.h
	ctags virus_gencor_simple.cc \
		poissonian.cc \
		poissonian.h

virus_gencor_simple: virus_gencor_simple.o poissonian.o
	g++ -o virus_gencor_simple virus_gencor_simple.o poissonian.o

virus_gencor_simple.o: virus_gencor_simple.cc
	g++ -c virus_gencor_simple.cc

virus_gencor: virus_gencor.o
	mpic++ -o virus_gencor virus_gencor.o

virus_gencor.o: virus_gencor.cc virus_gencor.h
	mpic++ -c virus_gencor.cc

virus_gencor_tauleap: virus_gencor_tauleap.o poissonian.o
	mpic++ -g -o virus_gencor_tauleap virus_gencor_tauleap.o poissonian.o

virus_gencor_tauleap.o: virus_gencor_tauleap.cc virus_gencor_tauleap.h poissonian.h
	mpic++ -g -c virus_gencor_tauleap.cc

poissonian.o: poissonian.cc poissonian.h
	mpic++ -g -c poissonian.cc

clean:
	rm -f virus_gencor \
		virus_gencor.o \
		virus_gencor_tauleap \
		virus_gencor_tauleap.o \
		poissonian.o
