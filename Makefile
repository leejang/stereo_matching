###############################
# Qulifying Exam Question 2
# Jangwon Lee
# leejang@indiana.edu
# November 24, 2016
###############################

# Compile and link flags
CC = g++

INCLUDES = -I./include
SRC_DIR = ./src
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

CFLAGS = -I./MRF_JW
LDFLAGS = -lpng -L./MRF_JW -lMRF

OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = stereo_reconstruction

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CC) $(INCLUDES) $(CFLAGS) -o $@ -c $< $(LDFLAGS)

clean:
	rm -rf $(SRC_DIR)/*.o $(EXECUTABLE)
