CC = mpicxx
ifeq ($(shell hostname), laura)
  CFLAGS = -I/usr/include -I/workspace/eigen-3.4.0 -L/usr/lib
else
  CFLAGS = -I/usr/include -I/usr/include/eigen3 -L/usr/lib
endif
#OFLAG = -fopenmp

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
