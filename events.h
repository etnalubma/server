#ifndef EVENTS_H
#define EVENTS_H

#include "enode.h"

typedef struct sevents * events;

/*
    Constructor de una lista dinamica de eventos.
*/
events create_events(void);

/*Destructor*/
events destroy_events(events s);

/*Agrega el evento `e` a la lista `s`*/
events insert_event(events s, event e);

/*Devuelve el evento mas proximo a ocurrir de la lista `s`*/
event next_event(events s);

void print_events(events s);

#endif
