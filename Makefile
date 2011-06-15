all: virus_gencor

virus_gencor: virus_gencor.o
	mpic++ -o virus_gencor virus_gencor.o

virus_gencor.o: virus_gencor.cc virus_gencor.h
	mpic++ -c virus_gencor.cc

virus_gencor_tauleap: virus_gencor_tauleap.o poissonian.o
	mpic++ -o virus_gencor_tauleap poissonian.o

virus_gencor_tauleap.o: virus_gencor_tauleap.cc poissonian.h
	mpic++ -c virus_gencor_tauleap.cc

poissonian.o: poissonian.cc poissonian.h
	mpic++ -c poissonian.cc

clean:
	rm -f virus_gencor \
		virus_gencor.o
