// poisson.h - Declares the function poisson() which returns a random
// number (cast as a double) belonging to a Poissonian distribution.

#ifndef POISSONIAN_H_
#define POISSONIAN_H_

double gammln(double);
double poissonian_knuth(double, unsigned short *);
double poissonian_reject(double, unsigned short *);
double poissonian(double, unsigned short *buf);

#endif /* POISSONIAN_H_ */
