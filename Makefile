CXX      = g++
CXXFLAGS =	-std=c++11 -O2 -DNDEBUG -g -Wall -fmessage-length=0
OBJS     =	$(addsuffix .o,$(addprefix src/,ms_index text_statistic bit-vector))

SDSL_INST= ~/arch/Darwin_x86_64
LIBS     = -L $(SDSL_INST)/lib -lsdsl -ldivsufsort -ldivsufsort64 
INCLUDE  = -I $(SDSL_INST)/include

TARGETS  = $(addprefix bin/,ms_index text_statistic bit-vector)

#$(TARGET):	$(OBJS)
#	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGETS)

$(TARGETS) : bin/% : src/%.o
	$(CXX) -o $@ $^ $(LIBS)

$(OBJS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)


clean:
	rm -f $(OBJS) $(TARGET)
