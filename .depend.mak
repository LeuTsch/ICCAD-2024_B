GlobalPlacer.o: GlobalPlacer.cpp GlobalPlacer.h solver.h inst.h parser.h
STAEngine.o: STAEngine.cpp inst.h parser.h STAEngine.h solver.h
inst.o: inst.cpp inst.h
legalizer.o: legalizer.cpp inst.h parser.h legalizer.h solver.h
main.o: main.cpp parser.h solver.h inst.h STAEngine.h legalizer.h
parser.o: parser.cpp parser.h
solver.o: solver.cpp inst.h solver.h parser.h STAEngine.h legalizer.h \
 GlobalPlacer.h
