#include <math.h>
#include <stdint.h>
#include <stdlib.h>

//This is an implementation of the exponential random number generation algorithm of Ahrens and Dieter (1972).
// It calls two functions generating a uniform random number:
// genrand_real3()  is expected to generate a random number in (0,1)
// genrand_real2()  is expected to generate a random number in [0,1)

//References:
//Computer methods for sampling from the exponential and normal distributions},
//J. H. Ahrens and U. Dieter, Communications of the ACM 15, 873 (1972)

//L. Devroye,
//Non-Uniform Random Variate Generation  (page 396 and following)
//Springer-Verlag, 1986
//http://www.nrbook.com/devroye/

void random_exponential_init()
{
    return;
}

double random_exponential()
{
    //First, we initialize some values
    //q[n] = \sum_{i=1}^{n+1} \log(2)^i / i!
    //(log is natural log)

    const static double q[] =
    {
       0.6931471805599453,
       0.9333736875190460,
       0.9888777961838676,
       0.9984959252914960,
       0.9998292811061389,
       0.9999833164100728,
       0.9999985691438767,
       0.9999998906925558,
       0.9999999924734159,
       0.9999999995283275,
       0.9999999999728814,
       0.9999999999985598,
       0.9999999999999290,
       0.9999999999999968,
       0.9999999999999999,
       1.
    };

    double u = genrand_real3();   //No zero (or 1) allowed!

    int exp;
    u = 2*frexp(u, &exp) - 1;

    if (u <= q[0]) return u - exp*q[0];  

    int i = 0;
    double ustar = genrand_real2(), umin = ustar;
    do {
        ustar = genrand_real2();
        if (umin > ustar)
            umin = ustar;
        i++;
    } while (u > q[i]);
    return (umin - exp) * q[0];
}

