#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "randomgen.h"
#include "estimate.h"
#include "server.h"
#include "sort.h"
#include "densities.h"
#include "tests.h"

#define BITS 32

float server_time(randgen rg);
double acum_gamma(double x);

int main(int argc, char *argv[]){
    randgen rg;
    server srv;
    FILE *file, *file2;
    int freturn = 1;
    float interval;
    int ocu;
    int *ocurrencies;
    double tattending, tmp = 0;
    int i, attended, sample;
    double t, a, s, p;
    int q, k;
    double fi, fj, t_zero, d_zero;
    double *results, *probabilities;   
    double mean, median, median1, median2, s2, skewness, min, max, ln_mean, ln_s2, alpha, beta, gamma_fn, deltab;
    
    /*Tiempo de atencion*/
    t = 8.;
    /*Razon de arrivos*/
    a = 4.;
    /*Razon de servicio*/
    s = 4.5;
    /*Cola de servidor*/
    q = 4;
    
    rg = create_rg(RAN2, 5, BITS);
    printf("\n");
    printf("q: %i, t: %f, a: %f, s: %f\n", q, t, a, s);                  
    printf("\n");    

    sample = 500;
    results = calloc(sample, sizeof(double));
    file = fopen("histogram.dat", "w");
    
    for(i=0;i<sample;i++){
        srv = create_server(rg, t, a, s, q);
        tmp = run_server(srv, &attended, &tattending);
        fprintf(file, "%f\n", tmp);
        results[i] = tmp;        
        
        srv = destroy_server(srv);
    }
    fclose(file);   
    
    file = fopen("scatter.dat", "w");
    for(i=0;i<sample-1;i++){
        fprintf(file, "%f %f\n", results[i], results[i+1]);
    }
    fclose(file);       


    printf("Estadisticos\n");
    printf("============\n");        
    bublesort(sample, results);

    min = results[0];
    max = results[sample-1];
    mean = estimate_mean(sample, results);
    ln_mean = estimate_ln_mean(sample, results);    
    median = estimate_median(sample, results);
    median1 = estimate_median(sample/2, results);
    median2 = estimate_median(sample/2, &results[sample/2]);
    s2 = estimate_s2(sample, results);
    ln_s2 = estimate_ln_s2(sample, results);    
    skewness = estimate_skewness(sample, results);
    
    printf("Minimo: %f\n", min);
    printf("Maximo: %f\n", max);    
    printf("Media: %f\n", mean);
    printf("Mediana: %f\n", median);
    printf("Mediana 1er cuartil: %f\n", median1);    
    printf("Mediana 2do cuartil: %f\n", median2);                
    printf("Varianza: %f\n", s2);
    printf("Asimetria: %f\n", skewness);
    printf("\n");
    
    file = fopen("boxplot.dat", "w");
    fprintf(file, "%i %f %f %f %f %f\n", 1, min, median1, median, median2, max);
    fclose(file);
    
    
    printf("Distribuciones\n");
    printf("==============\n");    
    printf("Normal(%f, %f)\n", mean, s2*((sample-1)/(double)sample));
    printf("Lognormal(%f, %f)\n", ln_mean, ln_s2);
    
    alpha = 10;/*9.539; /*Valor tabulado por T = 18.741333 = estimate_gamma_t(sample, results)*/

    gamma_fn = 362880; /*129979; /*integrate x^(alpha - 1)*e^(-x) from 0 to infinity*/
    beta = mean/alpha;
    printf("Gamma(%f, %f)\n", alpha, beta);    
    printf("\n");
    
    deltab = 0.0625;
    
    file = fopen("normal.dat", "w");    
    for(i=0;i<sample;i++){
        fprintf(file, "%f %f\n", results[i], deltab*normal(results[i], mean, s2));
    }
    fclose(file);  
    
    file = fopen("lognormal.dat", "w");    
    for(i=0;i<sample;i++){
        fprintf(file, "%f %f\n", results[i], deltab*lognormal(results[i], ln_mean, ln_s2));
    }
    fclose(file);  
    
    file = fopen("gamma.dat", "w");    
    for(i=0;i<sample;i++){
        fprintf(file, "%f %f\n", results[i], deltab*gamma(results[i], alpha, beta, gamma_fn));
    }
    fclose(file);  
    
    k = 25;
    probabilities = calloc(k, sizeof(double));
    ocurrencies = calloc(k, sizeof(int));
    
    file = fopen("histo.data", "r");    
    file2 = fopen("h.dat", "w");        
    for(i=0; i<k; i++){
        freturn = fscanf(file, "%f %i\n", &interval, &ocu);
        if(freturn!=EOF){
            p = ocu/(double)sample;
            fprintf(file2, "%f %f\n", interval, p);
        }
    }
    fclose(file);
    fclose(file2);
    /*
    printf("Normal\n");
    printf("======\n");

    file = fopen("histo.data", "r");    
    for(i=0; i<k-1; i++){
        freturn = fscanf(file, "%f %i\n", &interval, &ocu);
        if(freturn!=EOF){
            ocurrencies[i] = ocu;
            fi = cumulative_normal((i+1)*deltab, mean, s2*((sample-1)/(double)sample));
            fj = cumulative_normal(i*deltab, mean, s2*((sample-1)/(double)sample));
            probabilities[i] = fi - fj;
        }
    }
    fclose(file);
    
    t_zero = calculate_t(k, sample, ocurrencies, probabilities);
    printf("T cero %f\n", t_zero);
    printf("Chi2 %f\n", chi_square(k-1-2, t_zero));
    
    printf("LogNormal\n");
    printf("=========\n");
    
    file = fopen("histo.data", "r");    
    for(i=0; i<k-1; i++){
        freturn = fscanf(file, "%f %i\n", &interval, &ocu);
        if(freturn!=EOF){
            ocurrencies[i] = ocu;
            fi = cumulative_normal(log((i+1)*deltab), ln_mean, ln_s2);
            fj = cumulative_normal(log(i*deltab), ln_mean, ln_s2);
            probabilities[i] = fi - fj;
        }
    }
    fclose(file);    
    
    t_zero = calculate_t(k, sample, ocurrencies, probabilities);
    printf("T cero %f\n", t_zero);
    printf("Chi2 %f\n", chi_square(k-1-2, t_zero));    
*/
    printf("Gamma\n");
    printf("=====\n");
        
    file = fopen("histo.data", "r");
    file2 = fopen("gamma_acum.dat", "w");        
    for(i=0; i<k; i++){
        freturn = fscanf(file, "%f %i\n", &interval, &ocu);
        if(freturn!=EOF){
            ocurrencies[i] = ocu;
            fi = cumulative_gamma(i*deltab, alpha, beta);
            fj = cumulative_gamma((i-1)*deltab, alpha, beta);
            /*fprintf(file2, "%f %f\n", interval, fi-fj);*/
            probabilities[i] = fi - fj;
        }
    }
    fclose(file);    
    fclose(file2);    
    
    t_zero = calculate_t(k, sample, ocurrencies, probabilities);
    printf("T cero %f\n", t_zero);
    printf("Chi2 %f\n", chi_square(k-1-2, 66));    
    
    d_zero = calculate_d(sample, results, &acum_gamma);
    printf("D zero: %f\n", d_zero);
    printf("p-valor: %f\n", simulate_d(rg, 10000, sample, d_zero));
    free(probabilities);
    free(ocurrencies);
    free(results);         
        
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

double acum_gamma(double x){
    return cumulative_gamma(x, 10.000000, 0.050153);
}


