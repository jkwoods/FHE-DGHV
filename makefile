all: main.o Encoding.o Pk.o Deltas.o Pri_U.o utils.o PseudoRandomInts.o 
	g++ -o fhe main.o Encoding.o Pk.o Deltas.o Pri_U.o utils.o PseudoRandomInts.o -lgmpxx -lgmp -I/ccs/proj/gen119/woods/gmp-6.1.2 -L/ccs/proj/gen119/woods/gmp-6.1.2

main.o: main.cpp
	g++ -c main.cpp

Encoding.o: Encoding.cpp Encoding.hpp
	g++ -c Encoding.cpp

Pk.o: Pk.cpp Pk.hpp
	g++ -c Pk.cpp

Deltas.o: Deltas.cpp Deltas.hpp
	g++ -c Deltas.cpp

Pri_U.o: Pri_U.cpp Pri_U.hpp
	g++ -c Pri_U.cpp

utils.o: utils.cpp utils.hpp
	g++ -c utils.cpp

PseudoRandomInts.o: PseudoRandomInts.cpp PseudoRandomInts.hpp
	g++ -c PseudoRandomInts.cpp

clean:
	rm main.o Encoding.o Pk.o Deltas.o Pri_U.o utils.o PseudoRandomInts.o
