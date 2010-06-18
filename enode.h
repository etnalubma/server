#ifndef ENODE_H
#define ENODE_H

#include "event.h"

typedef struct senode * enode;

enode create_enode(event node);

enode destroy_enode(enode node);

enode get_next(enode node);

enode set_next(enode actual, enode next);

event get_event(enode node);

#endif
