# makefile template

adev : adev.cxx
	g++ -std=c++17 adev.cxx ../obj/toolbox.o -o ../bin/adev

../obj/toolbox.o : ../ToolBox/toolbox.cxx
	g++ -std=c++17 -c ../ToolBox/toolbox.cxx -o ../obj/toolbox.o
