# makefile template

TARGET=edev
	
edev : edev.cxx
	g++ -std=c++17 -g edev.cxx ../obj/toolbox.o -o ../bin/$(TARGET)
	
legendre : legendre.cxx
	g++ -std=c++17 -g legendre.cxx ../obj/toolbox.o -o ../bin/legendre
	
ddevloop : ddevloop.cxx
	g++ -std=c++17 -g ddevloop.cxx ../obj/toolbox.o -o ../bin/ddevloop
	
bdev : bdev.cxx
	g++ -std=c++17 -g bdev.cxx ../obj/toolbox.o -o ../bin/bdev
	
adev : adev.cxx
	g++ -std=c++17 -g adev.cxx ../obj/toolbox.o -o ../bin/adev

../obj/toolbox.o : ../ToolBox/toolbox.cxx
	g++ -std=c++17 -c ../ToolBox/toolbox.cxx -o ../obj/toolbox.o

phoney : run gdb

run:
	../bin/$(TARGET)
	
gdb:
	gdb ../bin/$(TARGET)
