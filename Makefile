default: build run

build:
	gcc -lm -lGLU -lGL -lX11 c3d.c tigr.c -Ofast -flto -fopenmp -mavx -o c3d

run:
	./c3d models/bust.obj textures/wood.ppm
