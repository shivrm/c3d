build:
	gcc -lm -lGLU -lGL -lX11 c3d.c tigr.c -O2 -o c3d
	./c3d models/teapot.obj
