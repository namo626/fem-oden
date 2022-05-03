CC = g++
CFLAGS = -DEIGEN_DONT_PARALLELIZE -I/usr/include -I/usr/include/eigen3 -L/usr/lib

main: main.cpp Element.o Global.o Post.o
	$(CC) -o $@ $^ $(CFLAGS) $(OFLAG)

Element.o: Element.cpp
	$(CC) -c $< $(CFLAGS) $(OFLAG)

Global.o: Global.cpp
	$(CC) -c $< $(CFLAGS)

Post.o: Post.cpp
	$(CC) -c $< $(CFLAGS)

clean:
	rm -rf *.o main
