MFEM_DIR = $(HOME)/Documents/software/mfem/mfem_package/mfem-4.3

include $(MFEM_DIR)/config/config.mk

#gnu compiler
#CXX		= $(MFEM_CXX)
#intel compiler
CXX		= mpicxx

#Ubuntu
CXXFLAGS	= $(MFEM_FLAGS) -no-pie

#Centos
#CXXFLAGS	= $(MFEM_FLAGS)

INCLUDES	= $(wildcard src/*.h)
SRCS		= $(wildcard src/*.cpp)
OBJS		= $(patsubst %.cpp, %.o, $(SRCS))

Ocean_INCLUDE	= -I src

%.o : %.cpp
	@echo "GlobalMT is compiling "$<"..."
	@$(CXX) $(CXXFLAGS) $(Ocean_INCLUDE) -c $< -o $@ 

GlobalMT: $(OBJS)
	@$(CXX) $(CXXFLAGS) -o GlobalMT $(OBJS) $(MFEM_LIBS)

clean:
	@rm -rf $(OBJS)
	
cleanall: 
	@rm -rf $(OBJS) GlobalMT


# mpirun -np 1 ../../GlobalMT M2_tide.config |tee a.log
