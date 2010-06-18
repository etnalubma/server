#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "events.h"

struct sevents{
    enode first;
};

events create_events(void){
    events s = (events)calloc(1, sizeof(struct sevents));
    s->first = NULL;
    return s;
}

events destroy_events(events s){
    enode tmp;
    while(s->first){
        tmp = s->first;
        s->first = get_next(tmp);
        
        destroy_event(get_event(tmp));
        destroy_enode(tmp);
    }
    free(s);
    s = NULL;
    return s;
}

events insert_event(events s, event e){
    enode new_node, tmp_node, prev_node;
    event tmp_event;
    double t, new_time;
    
    assert(s && e);

    new_node = create_enode(e);
    new_time = get_time(e);
    tmp_node = s->first;
    
    if(!tmp_node){
        /*si s esta vacio*/
        s->first = new_node;
    }
    else{
        /*si s no esta vacio*/
        tmp_event = get_event(tmp_node);
        t = get_time(tmp_event);
        if( new_time <= t){
            /*se compara con el primer nodo*/
            new_node = set_next(new_node, tmp_node);
            s->first = new_node;            
        }
        else{
            /*se busca la posicion indicada*/
            prev_node = tmp_node;
            tmp_node = get_next(prev_node);
            while(tmp_node){
                tmp_event = get_event(tmp_node);
                t = get_time(tmp_event);
                if(new_time <= t)
                    break;
                else{
                    prev_node = tmp_node;
                    tmp_node = get_next(prev_node);
                }
            }
            set_next(prev_node, new_node);
            set_next(new_node, tmp_node);
        }        
    }
    
    return s;
}

event next_event(events s){
    enode tmp_node;
    event tmp_event;
    event e = NULL;
    double t;

    if(s->first){
        tmp_node = s->first;
        e = get_event(tmp_node);
        t = get_time(e);
       
        s->first = get_next(tmp_node);
        
        destroy_enode(tmp_node);

        tmp_node = s->first;
        while(tmp_node){
            /*se recorren todos los nodos para actualizar los tiempos*/
            tmp_event = get_event(tmp_node);
            tmp_event = update_time(tmp_event, t);
            
            tmp_node = get_next(tmp_node);
        }
    }
    
    return e;
}

