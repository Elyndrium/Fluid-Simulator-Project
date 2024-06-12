
build:
	del sim.exe
	g++ main.cpp -Ofast -flto -funroll-loops -finline-functions -march=native -o sim.exe

run: clean build
	sim.exe

debug:
	del debug_sim.exe
	g++ main.cpp -Og -Wall -Wextra -Wpedantic -I .stb_image_write.h -o debug_sim.exe

clean:
	del sim.exe
	del debug_sim.exe
	del *.png

tutte:
	del tutte.exe
	g++ tutte_mapping.cpp -O3 -o tutte.exe