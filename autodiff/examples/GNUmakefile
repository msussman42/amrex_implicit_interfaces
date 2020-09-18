CXX	=g++
CXXFLAGS=-g

# Normal object files to be created
OBJECTS = test_AD.o
# Object files created with no automatic differentiation
#NO_AD_OBJECTS = test.o
# Program name
PROGRAM = test_AD
# Include-file location
INCLUDES = -I/usr/local/include -I../adept-2.0.5/include
# Library location and name, plus the math library
#LIBS = -L/usr/local/lib -lm -ladept
LIBS = -L/usr/local/lib -lm 
# Rule to build the program (typing "make" will use this rule)
$(PROGRAM): $(OBJECTS) $(NO_AD_OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(NO_AD_OBJECTS) $(LIBS) -o $(PROGRAM)
# Rule to build a normal object file (used to compile all objects in OBJECTS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<
# Rule to build a no-automatic-differentiation object (used to compile ones in NO_AD_OBJECTS)
%_noad.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DADEPT_NO_AUTOMATIC_DIFFERENTIATION -c $< -o $@
