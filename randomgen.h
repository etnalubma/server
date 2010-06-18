#ifndef RANDOMGEN_H
#define RANDOMGEN_H

typedef enum GEN_METHOD{
    RAN2, 
    MZRAN13
} genmethod;

typedef struct srandgen *randgen;

/*
    Constructor del generador de numeros aleatorios.
    
    Parametros:
    @method RAN2 | MZRAN13
    @seed: Semilla para RAN2 (no se usa con MZRAN13)
    @bits: 32 | 64 (no se usa con RAN2)
*/
randgen create_rg(genmethod method, long seed, int bits);

/*Destructor*/
randgen destroy_rg(randgen rg);

/*Genera un numero aleatorio*/
float random_gen(randgen rg);

#endif
