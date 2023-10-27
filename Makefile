all: lab3

CC=oshcc
CFLAGS=-Wno-unused-parameter -Wno-format-overflow -Wall -Wextra -Werror -g -rdynamic -O3 -msse4.2
#CFLAGS=-fsanitize=address -fsanitize=undefined -Wno-unused-parameter -Wno-format-overflow -Wall -Wextra -Werror -g -rdynamic -O3 -msse4.2

lab3: lab3.o wctimer.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

nbody: 
	clang -o nbody nbody.c -lm

clean:
	rm -f *.o lab3 nbody
