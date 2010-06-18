#ifndef SERVER_H
#define SERVER_H

#include "randomgen.h"

typedef struct sserver * server;

/*
    @rg: Generador de v.a.
    @t: Horas en servicio por dia.
    @c: Razon de llegada.
    @s: Razon de servicio.
    @n: Tam. de la cola de espera.
*/
server create_server(randgen rg, double t, double c, double s, int n);

server destroy_server(server s);

double run_server(server srv, int *attended, double *tattending);

#endif

