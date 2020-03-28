X = g++
CXXFLAGS = -Wall -I/usr/include/eigen3/

proj1: Scene.o main.cpp
	$(CXX) $(CXXFLAGS) -O3 Scene.o main.cpp -o proj3
  
Scene.o: Scene.cpp Scene.h
	$(CXX) $(CXXFLAGS) -c Scene.cpp

clean:
	rm *.o*
	rm *~ 

run:
	./proj3
