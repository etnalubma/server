#include <assert.h>
#include <stdlib.h>

#include "enode.h"

struct senode{
    event e;
    enode next;
};

enode create_enode(event e){    
    enode node;

    assert(e);    
    
    node = (enode)calloc(1, sizeof(struct senode));
    node->e = e;
    node->next = NULL;
    return node;
}

enode destroy_enode(enode node){
    free(node);
    node = NULL;
    return node;
}

enode get_next(enode node){
    assert(node);
    return node->next;
}

enode set_next(enode actual, enode next){
    assert(actual);
    actual->next = next;
    return actual;
}

event get_event(enode node){
    assert(node);
    return node->e;
}

