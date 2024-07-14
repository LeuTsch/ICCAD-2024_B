CSRCS     = $(wildcard *.cpp)
CHDRS     = $(wildcard *.h)
COBJS     = $(addsuffix .o, $(basename $(CSRCS)))

CXX       = g++
CFLAGS    = -std=c++11 -g -O3 -Wall
EXTINCDIR = 
LIBDIR    = 
ECHO      = /bin/echo

EXEC      = main

.PHONY: depend

$(EXEC): $(COBJS)
	@$(ECHO) "> building: $@"
	@$(CXX) -o $@ $(CFLAGS) $(COBJS)

%.o : %.cpp
	@$(ECHO) "> compiling: $<"
	@$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	@rm -f $(COBJS) $(EXEC)

depend: .depend.mak
.depend.mak: $(CSRCS) $(CHDRS)
	@$(ECHO) Making dependencies ...
	@$(CXX) -MM $(DEPENDDIR) $(CSRCS) > $@

include .depend.mak