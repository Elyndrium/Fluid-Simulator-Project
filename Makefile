
build:
	del sim.exe
	g++ main.cpp -O3 -o sim.exe

# -Ofast -flto -funroll-loops -finline-functions -march=native

run: clean build
	sim.exe

debug:
	del debug_sim.exe
	g++ main.cpp -Og -Wall -Wextra -Wpedantic -I .stb_image_write.h -o debug_sim.exe

clean:
	del sim.exe
	del debug_sim.exe
	del *.svg
	del *.png
