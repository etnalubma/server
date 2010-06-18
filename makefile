CC = gcc
CFLAGS = -Wall -ansi -pedantic -g

OBJS = randomgen.o generators.o estimate.o event.o enode.o events.o server.o main.o

main: ${OBJS}
	${CC} -o program -lm ${OBJS}
		
clean:
	rm -f ${OBJS} program *~ *.dat
        
