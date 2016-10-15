CXX      = g++
CXXFLAGS = -std=c++11 -O2 -DNDEBUG -g -Wall -fmessage-length=0
OBJS     =	$(addsuffix .o,$(addprefix bin/,ms_index text_statistic bit-vector \
										range-minimum-queries fabio_djamal_ms))

SDSL_INST= ~/arch/Darwin_x86_64
LIBS     = -L $(SDSL_INST)/lib -lsdsl -ldivsufsort -ldivsufsort64 
INCLUDE  = -I $(SDSL_INST)/include

OBJDIR   = objects
TARGETS  = $(addprefix bin/,ms_index text_statistic bit-vector range-minimum-queries\
				fabio_djamal_ms)


all:	$(TARGETS)

$(TARGETS) : bin/% : bin/%.o
	$(CXX) -o $@ $^ $(LIBS)

$(OBJS) : bin/%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)


clean:
	rm -f $(OBJS) $(TARGET)
