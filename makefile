all: main.o Encoding.o Pk.o utils.o PseudoRandomInts.o 
	g++ -o fhe main.o Encoding.o Pk.o utils.o PseudoRandomInts.o -lgmpxx -lgmp -fopenmp -I/ccs/proj/gen119/woods/gmp-6.1.2 -L/ccs/proj/gen119/woods/gmp-6.1.2

main.o: main.cpp
	g++ -c main.cpp -fopenmp

Encoding.o: Encoding.cpp Encoding.hpp
	g++ -c Encoding.cpp -fopenmp

Pk.o: Pk.cpp Pk.hpp
	g++ -c Pk.cpp -fopenmp

utils.o: utils.cpp utils.hpp
	g++ -c utils.cpp -fopenmp

PseudoRandomInts.o: PseudoRandomInts.cpp PseudoRandomInts.hpp
	g++ -c PseudoRandomInts.cpp -fopenmp

clean:
	rm main.o Encoding.o Pk.o utils.o PseudoRandomInts.o
