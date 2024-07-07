STAEngine.o: STAEngine.cpp inst.h parser.h solver.h STAEngine.h
inst.o: inst.cpp inst.h
main.o: main.cpp parser.h solver.h inst.h STAEngine.h
parser.o: parser.cpp parser.h
solver.o: solver.cpp inst.h solver.h parser.h
