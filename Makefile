#   
#	File name		: Makefile
#	Date			: 
#	Version			: 
#	Author			: 
#

DEST			= .

HDRS			= global.hpp					\
				  src/initialization/initialization.hpp		\
				  src/initializationfunctions.hpp				\
				  src/neighbor/neighbor.hpp				\
				  src/diffusion/dc_operator.hpp			\
				  src/advection/advection.hpp				\
				  src/remeshing/remeshing.hpp				\

LIBS			= -larmadillo

INPS			= 

COMPILER		= g++ 

OPTFLAG			= -std=c++11 -O2
 
MAKEFILE		= Makefile


PROGRAM			= Test_Vortex

SRCS			= main.cpp						\
				  src/initialization/init_uniform.cpp		\
  				  src/initialization/init_random.cpp		\
				  src/initialization/create_index.cpp		\
				  src/initialization/functions.cpp				\
				  src/neighbor/create_neighbor_list.cpp				\
				  src/neighbor/create_neighbor_list_inter.cpp				\
				  src/neighbor/ghost_particle.cpp				\
				  src/neighbor/ghost_particle_unsym.cpp				\
				  src/neighbor/ghost_particle_unsym_v2.cpp				\
				  src/neighbor/neighbor_between.cpp				\
				  src/neighbor/neighbor_temp.cpp				\
				  src/diffusion/dc_pse.cpp			\
  				  src/diffusion/vandermonde.cpp			\
  				  src/diffusion/vandermonde_zero.cpp			\
				  src/advection/advection.cpp				\
				  src/remeshing/dc_remesh.cpp				\
				  src/remeshing/m4_remesh.cpp				\

OBJS			= $(SRCS:.cpp=.o) 	

.cpp.o:
			$(COMPILER) $(OPTFLAG) -c $*.cpp -o $*.o 

all:			$(PROGRAM)

$(PROGRAM):		$(OBJS) $(LIBS)
				@echo -n "Loading Program $(PROGRAM) ... "
				@$(COMPILER) $(OPTFLAG) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
				@echo "done"

clean:;			@rm -f $(SRCS:.cpp=.o) $(SRCS:.cpp=.il) $(PROGRAM)



