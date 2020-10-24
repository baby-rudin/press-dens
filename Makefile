press-dens: *.o
	g++ -O3 -o press-dens *.o

*.o: *.cpp
	g++ -O3 -c *.cpp -o *.o

.PHONY: clean
clean:
	rm -rf *.o press-dens