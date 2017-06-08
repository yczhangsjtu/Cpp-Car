car: car.cpp drawer.h
	g++ $< -o $@ -lGL -lGLU -lglut -l3ds
