#include <math.h>
#include "generators.h"

double exponential(randgen rg, double lambda){
    float rand;

    rand = random_gen(rg);
    return -(1/lambda)*log(rand);
}
