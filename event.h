#ifndef EVENT_H
#define EVENT_H

typedef enum EVENT_TYPE {
    IN,
    OUT
} event_type;

typedef struct sevent * event;

/*
    Constructor de eventos.
    Parametros:
    @time Tiempo que falta para que ocurra el evento
    @type IN | OUT
*/
event create_event(double time, event_type type);

/*Destructor*/
event destroy_event(event e);

event_type get_type(event e);

double get_time(event e);

/*Actualiza el tiempo del evento*/
event update_time(event e, double time);

#endif
