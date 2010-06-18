#include <assert.h>
#include <stdlib.h>

#include "event.h"

struct sevent{
    double time;
    event_type type;
};

event create_event(double time, event_type type){
    event e;
    
    e = (event)calloc(1, sizeof(struct sevent));
    e->time = time;
    e->type = type;
    return e;
}

event destroy_event(event e){
    free(e);
    e = NULL;
    return e;
}

event_type get_type(event e){
    assert(e);
    return e->type;
}

double get_time(event e){
    assert(e);
    return e->time;
}

event update_time(event e, double time){
    assert(e);
    assert(e->time > time);
    
    e->time = e->time - time;
    return e;
}

