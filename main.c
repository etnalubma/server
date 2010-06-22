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
#include "special.h"
#define BITS 32

void run_server_simulation(randgen rg, double *results, int n);
void write_histogram(double *results, int n);
void write_scatter(double *results, int n);
void write_boxplot(double min, double max, double median, double median1, double median2);

void normal_frequencies(double *results, int n, double mean, double s2, double deltab);
void lognormal_frequencies(double *results, int n, double mean, double s2, double deltab);
void gamma_frequencies(double *results, int n, double alpha, double beta, double gamma_fn, double deltab);

void calculate_probabilities_gamma(int k, double deltab, double alpha, double beta, double gamma_fn, double *probabilities);
void calculate_probabilities_lognormal(int k, double deltab, double ln_mean, double ln_s2, double *probabilities);
void calculate_probabilities_normal(int k, double deltab, double mean, double s2, double *probabilities);
void calculate_ocurrencies(double *results, int n, int k, double deltab, int *ocurrencies);

double cumulative_gamma(double x);
double cumulative_lognormal(double x);
double cumulative_normal(double x);

int main(int argc, char *argv[]){
    randgen rg;

    int *ocurrencies;

    int sample;
    int k;
    double t_zero, d_zero;
    double *results, *probabilities;   
    double mean, median, median1, median2, s2, skewness, min, max, ln_mean, ln_s2, alpha, beta, gamma_fn, deltab;
        
    rg = create_rg(RAN2, 1, BITS);
    printf("\n");        

    sample = 500;
    results = calloc(sample, sizeof(double));
    
    run_server_simulation(rg, results, sample);

    /*Actividad 1*/    
    write_scatter(results, sample);
    
    
    /*Actividad 2*/
    printf("Estadisticos\n");
    printf("============\n");        
    
    /*2.b)*/
    write_histogram(results, sample);   
        
    bublesort(sample, results);

    /*2.a)*/
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

    /*2.c)*/
    write_boxplot(min, max, median, median1, median2);
            
    /*Actividad 3*/
        
    printf("Distribuciones\n");
    printf("==============\n");    
    printf("Normal(%f, %f)\n", mean, s2*((sample-1)/(double)sample));
    printf("Lognormal(%f, %f)\n", ln_mean, ln_s2);
    
    alpha = 9.789; /*Valor tabulado por T = 19.110559 = estimate_gamma_t(sample, results)*/
    gamma_fn = 226175; /*integrate x^(alpha - 1)*e^(-x) from 0 to infinity*/
    beta = mean/alpha;
    printf("Gamma(%f, %f)\n", alpha, beta);    

    printf("\n");    
    
    /*Actividad 4*/
    k = 25;
    deltab = (max+0.01)/k;
    
    /*4.a)*/    
       
    normal_frequencies(results, sample, mean, s2*((sample-1)/(double)sample), deltab);
    lognormal_frequencies(results, sample, ln_mean, ln_s2, deltab);
    gamma_frequencies(results, sample, alpha, beta, gamma_fn, deltab);
    
    /*4.b)*/
    
    ocurrencies = calloc(k, sizeof(int));    
    calculate_ocurrencies(results, sample, k, deltab, ocurrencies);
    
    printf("\nNormal\n");
    printf("======\n");
    
    probabilities = calloc(k, sizeof(double));    

    calculate_probabilities_normal(k, deltab, mean, s2*((sample-1)/(double)sample), probabilities);    

    t_zero = calculate_t(k, sample, ocurrencies, probabilities);    
    printf("T cero %f\n", t_zero);
    printf("Chi2 %f\n", chi_square(k-1-2, t_zero));    
    
    d_zero = calculate_d(sample, results, &cumulative_normal);
    printf("D zero: %f\n", d_zero);
    printf("p-valor: %f\n", simulate_d(rg, 1000, sample, d_zero));
    
    free(probabilities);    
        
    printf("\nLog Normal\n");
    printf("==========\n");
    probabilities = calloc(k, sizeof(double));    

    calculate_probabilities_lognormal(k, deltab, ln_mean, ln_s2, probabilities);

    t_zero = calculate_t(k, sample, ocurrencies, probabilities);
    printf("T cero %f\n", t_zero);
    printf("Chi2 %f\n", chi_square(k-1-2, t_zero));
       
    d_zero = calculate_d(sample, results, &cumulative_lognormal);
    printf("D zero: %f\n", d_zero);
    printf("p-valor: %f\n", simulate_d(rg, 1000, sample, d_zero));
    
    free(probabilities);    

    printf("\nGamma\n");
    printf("=====\n");
    probabilities = calloc(k, sizeof(double));
    
    calculate_probabilities_gamma(k, deltab, alpha, beta, gamma_fn, probabilities);
    
    t_zero = calculate_t(k, sample, ocurrencies, probabilities);    
    printf("T cero %f\n", t_zero);
    printf("Chi2 %f\n", chi_square(k-1-2, t_zero));    
    
    d_zero = calculate_d(sample, results, &cumulative_gamma);
    printf("D zero: %f\n", d_zero);
    printf("p-valor: %f\n", simulate_d(rg, 1000, sample, d_zero));
  
    free(probabilities);
    
    free(ocurrencies);
    free(results);             
    rg = destroy_rg(rg);
    printf("\n");
    return 0;
}

void run_server_simulation(randgen rg, double *results, int n){
    double t, a, s;
    double tattending;
    int attended;
    int q, i;
    
    server srv;
    
    /*Tiempo de atencion*/
    t = 8.;
    /*Razon de arrivos*/
    a = 4.;
    /*Razon de servicio*/
    s = 4.5;
    /*Cola de servidor*/
    q = 4;
    
    printf("\n");
    printf("q: %i, t: %f, a: %f, s: %f\n", q, t, a, s);                  
    printf("\n");    
    
    for(i=0;i<n;i++){
        srv = create_server(rg, t, a, s, q);
        results[i] = run_server(srv, &attended, &tattending);
        srv = destroy_server(srv);        
    }
}

void write_histogram(double *results, int n){
    FILE *file;
    int i;
    file = fopen("histogram.dat", "w");
    
    for(i=0;i<n;i++){
        fprintf(file, "%f\n", results[i]);
    }
    
    fclose(file);
}   

void write_scatter(double *results, int n){
    FILE *file;
    int i;
    file = fopen("scatter.dat", "w");
    
    for(i=0;i<n-1;i++){
        fprintf(file, "%f %f\n", results[i], results[i+1]);
    }
    fclose(file);
}

void write_boxplot(double min, double max, double median, double median1, double median2){
    FILE *file;
    
    file = fopen("boxplot.dat", "w");
    fprintf(file, "%i %f %f %f %f %f\n", 1, min, median1, median, median2, max);
    fclose(file);
}

void normal_frequencies(double *results, int n, double mean, double s2, double deltab){
    FILE *file;
    int i;
    file = fopen("normal.dat", "w");    
    for(i=0;i<n;i++){
        fprintf(file, "%f %f\n", results[i], n*deltab*normal(results[i], mean, s2));
    }
    fclose(file);
}

void lognormal_frequencies(double *results, int n, double mean, double s2, double deltab){
    FILE *file;
    int i;
    file = fopen("lognormal.dat", "w");    
    for(i=0;i<n;i++){
        fprintf(file, "%f %f\n", results[i], n*deltab*lognormal(results[i], mean, s2));
    }
    fclose(file);  
}

void gamma_frequencies(double *results, int n, double alpha, double beta, double gamma_fn, double deltab){
    FILE *file;    
    int i;
    file = fopen("gamma.dat", "w");    
    for(i=0;i<n;i++){
        fprintf(file, "%f %f\n", results[i], n*deltab*gamma(results[i], alpha, beta, gamma_fn));
    }
    fclose(file);  
}

void calculate_ocurrencies(double *results, int n, int k, double deltab, int *ocurrencies){
    int i, j;
    for(i=0; i<k; i++){
        for(j=0;j<n;j++) {
            if (results[j] > (deltab*i) && results[j] <= (deltab*(i+1)))
                ocurrencies[i]++;
        }
    }
}

void calculate_probabilities_normal(int k, double deltab, double mean, double s2, double *probabilities){
    int i;
    double m, y;
    m = deltab/2.;
    for(i=0; i<k; i++){
        y = i*deltab + m;
        probabilities[i] = deltab*normal(y, mean, s2);
    }

}

void calculate_probabilities_lognormal(int k, double deltab, double ln_mean, double ln_s2, double *probabilities){
    int i;
    double m, y;
    m = deltab/2.;
    for(i=0; i<k; i++){
        y = i*deltab + m;
        probabilities[i] = deltab*lognormal(y, ln_mean, ln_s2);
    }

}

void calculate_probabilities_gamma(int k, double deltab, double alpha, double beta, double gamma_fn, double *probabilities){
    int i;
    double m, y;
    m = deltab/2.;
    for(i=0; i<k; i++){
        y =  i*deltab + m;
        probabilities[i] = deltab*gamma(y, alpha, beta, gamma_fn);
    }
}

double cumulative_normal(double x){
    double y, mean, s2;

    mean = 0.509307;
    s2 = 0.026826;
    
    y = (x-mean)/sqrt(s2);
    return 0.5*(1 + erff(y/sqrt(2.)));
}

double cumulative_lognormal(double x){
    double ln_mean, ln_s2, y;

    ln_mean = -0.727032;
    ln_s2 = 0.108136;

    y = (log(x) - ln_mean)/sqrt(ln_s2); /*log(x)~N(ln_mean, ln_s2) => y ~ N(0,1)*/
    return 0.5*(1 + erff(y/sqrt(2.)));
}

double cumulative_gamma(double x){
    double alpha, beta, gamma_fn;
    
    alpha = 9.789000;
    beta = 0.052028;
    gamma_fn = 226175;
    
    return gammp(alpha, x/beta);
}
