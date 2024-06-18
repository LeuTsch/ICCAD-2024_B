inst.o: inst.cpp inst.h
main.o: main.cpp parser.h
parser.o: parser.cpp parser.h
solver.o: solver.cpp inst.h solver.h parser.h
