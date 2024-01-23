# makefile template

bdev : bdev.cxx
	g++ -std=c++17 -g bdev.cxx ../obj/toolbox.o -o ../bin/bdev
	
adev : adev.cxx
	g++ -std=c++17 -g adev.cxx ../obj/toolbox.o -o ../bin/adev

../obj/toolbox.o : ../ToolBox/toolbox.cxx
	g++ -std=c++17 -c ../ToolBox/toolbox.cxx -o ../obj/toolbox.o

phoney : run

run:
	../bin/bdev
