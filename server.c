#include <stdlib.h>
#include <assert.h>

#include "events.h"
#include "generators.h"
#include "server.h"


struct sserver{
    events events_list;
    double t;
    double c;
    double s;
    int n;
    randgen rg;
};

server create_arrivals(server srv);

server create_server(randgen rg, double t, double c, double s, int n){
    server srv;
    
    srv = (server)calloc(1, sizeof(struct sserver));
    
    srv->events_list = create_events();
    srv->t = t;
    srv->c = c;
    srv->s = s;
    srv->n = n;
    srv->rg = rg;
    
    return srv;
}

server destroy_server(server srv){
    srv->events_list = destroy_events(srv->events_list);
    free(srv);
    srv = NULL;
    return srv;
}

double run_server(server srv, int *attended, double *tattending){
    event current, e;
    
    double t, tatt=0, wait=0, trunning=0;
    int att=0, q = 0;
    
    /*Crea los tiempos de arrivo al servidor*/
    srv = create_arrivals(srv);

    while(trunning < srv->t){ 
        /*Hasta que se alcance el tiempo de atencion del servidor*/
        current = next_event(srv->events_list);
        trunning += get_time(current);
        
        if(q > 0){
            /*Si hay clientes en la cola del servidor, se decrementa el tiempo
            de espera con el tiempo transcurrido hasta el evento `current`
            */
            wait -= get_time(current);
        }
                        
        if(get_type(current) == IN){
            /*El evento `current` es un cliente que quiere ingresar al servidor*/
            if(q < srv->n){
                /*Si hay lugar en la cola del servidor, se agrega el cliente.
                De lo contrario, se descarta.
                */
                q++;
                t = exponential(srv->rg, srv->s);
                
                /*Se actualiza el tiempo de espera en la cola*/
                wait += t;
                
                /*Se actualiza el tiempo total de clientes en la cola 
                y la cantidad de clientes atendidos.
                */
                tatt += wait;
                att++;
                
                e = create_event(wait, OUT);
                srv->events_list = insert_event(srv->events_list, e);
            }
        }
        else{
            /*El evento `current` es un cliente que ya fue atendido por el servidor,
            luego se libera un lugar en la cola.
            */
            q--;     
        }
        destroy_event(current);
    }
    *attended = att;
    *tattending = tatt;
    return tatt/att;
}

server create_arrivals(server srv){
    /*Se generan los tiempos de arrivo hasta alcanzar el tiempo
    maximo de atencion.
    */
    event e;
    double t=0;
    
    while(t < srv->t){
        t += exponential(srv->rg, srv->c);
        e = create_event(t, IN);
        srv->events_list = insert_event(srv->events_list, e);
    }
    
    return srv;
}
