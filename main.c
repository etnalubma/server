#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "randomgen.h"
#include "estimate.h"
#include "server.h"

#define BITS 32

int confidence_interval95(float d, float s, int i);
float server_time(randgen rg);

int uni_gen_interval(randgen rg, int from, int to){
    float rand = random_gen(rg);
    return (int)(rand*(to-from+1)) + from;
}


int main(int argc, char *argv[]){
    randgen rg;
    server srv;
    estimate est;
    int z, sample;
    double theta, d, n, tattending, tmp = 0, D, N;
    int iterations=100000, i, j, k, r=100, attended;
    double t, a, s, y;
    int q;
    double *times, *clients;
    
    float mean=0, sigma=0;
    int samples=0;
    
    /*Tiempo de atencion*/
    t = 8.;
    /*Razon de arrivos*/
    a = 4.;
    /*Razon de servicio*/
    s = 4.2;
    /*Cola de servidor*/
    q = 4;
    
    rg = create_rg(RAN2, 5, BITS);
    printf("\n");
    printf("q: %i, t: %f, a: %f, s: %f\n", q, t, a, s);                  
    printf("\n");    


    d = 0;
    n = 0;
    
    sample = 100;
    times = calloc(sample, sizeof(double));
    clients = calloc(sample, sizeof(double));
    
    for(i=0;i<sample;i++){
        srv = create_server(rg, t, a, s, q);
        tmp = run_server(srv, &attended, &tattending);
        
        times[i] = tattending;
        clients[i] = attended;
        
        d += tattending;
        n += attended;
        
        srv = destroy_server(srv);
    }
    theta = d/n;
    
    printf("%i dias: theta = %f\n\n", sample, theta);    
    printf("NS\t ECM\n");              
    
    for(j=1; j<iterations; j*=10){
        y = 0;
        for(k=0; k<r; k++){
            D=0;
            N=0;
            for(i=0; i<j; i++){
                z = uni_gen_interval(rg, 0, sample-1);
                D += times[z];
                N += clients[z];
            }
            y += pow((D/N) - theta, 2)/(double)r;
        }       
        printf("\n");
        printf("%i \t %f\n", j, y);
    }
    free(times);
    free(clients);
    
    printf("\n");
    printf("======== Intervalo de confianza del 95%% y ancho de intervalo 0.04 ========\n");           
    printf("\n");
    est = estimateInit(rg, &server_time, 30, 0, &confidence_interval95);
    
    sample_mean_sigma(est, &mean, &sigma, &samples);
    
    printf("Media muestral: %f\n", mean);
    printf("Muestras (dias simulados): %i\n", samples);
    
    est = estimateDestroy(est);
    printf("\n");
    
    rg = destroy_rg(rg);

    return 0;
}

float server_time(randgen rg){
    server srv;
    double t, tattending;
    int attended;
    
    srv = create_server(rg, 8., 4., 4.2, 4);
    t = run_server(srv, &attended, &tattending);
    srv = destroy_server(srv);
    return t;
}

int confidence_interval95(float d, float s, int i) {
    return confidence_interval(s/sqrt(i), Z95, 0.04);
}

