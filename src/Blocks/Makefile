# Makefile - Blocks Folder

TARGET=classblock

$(TARGET): $(TARGET).cxx
	g++ -std=c++17 -pthread -g $(TARGET).cxx ../../obj/toolbox.o -o ../../bin/$(TARGET)
	
#-------------------------------------------------------------------------------

logbase: logbase.cxx
	g++ -std=c++17 -g logbase.cxx ../../obj/toolbox.o -o ../../bin/logbase
	
s3000: S3000.cxx
	g++ -std=c++17 -g S3000.cxx ../../obj/toolbox.o -o ../../bin/modS3000
	
#-------------------------------------------------------------------------------
phoney : run gdb

run:
	../../bin/$(TARGET)
	
gdb:
	gdb ../../bin/$(TARGET)
