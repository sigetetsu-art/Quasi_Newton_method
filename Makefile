CC = g++
CFLAG = -O4 -Wall

all: quasi_newton

quasi_newton: Quasi_newton.cpp
						$(CC) $(CFLAG) Quasi_newton.cpp -o quasi_newton

clean:
		rm -f quasi_newton *.o*~