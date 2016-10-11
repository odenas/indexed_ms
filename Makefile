CXX      = g++
CXXFLAGS =	-std=c++11 -O2 -DNDEBUG -g -Wall -fmessage-length=0
OBJS     =	src/ms_index.o

SDSL_INST= ~/arch/Darwin_x86_64
LIBS     = -L $(SDSL_INST)/lib -lsdsl -ldivsufsort -ldivsufsort64 
INCLUDE  = -I $(SDSL_INST)/include
TARGET   = bin/ms_index

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
