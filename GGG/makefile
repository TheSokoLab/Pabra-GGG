#------------------------------#
# NEW MAKEFILE                 #
#                              #
# made by Tomek, Mar 2010      #
#------------------------------#

# Some variables used in this file, set to empty by default
# The values are set later for every compilation rule specifically
NT =
OPT64 =
REMEMBER =

# Source and object files
TOOLS = ./Tools/
HF = $(TOOLS)mersenne_twister.hpp $(TOOLS)array2D.hpp $(TOOLS)tree_node.hpp $(TOOLS)binary_tree.hpp $(TOOLS)random_variable.hpp $(TOOLS)histo1D.hpp $(TOOLS)histo2D.hpp run_parameter.hpp parameters.hpp reaction.hpp observable.hpp nucleus.hpp cortex.hpp simulation.hpp
OF = random_variable.o histo1D.o histo2D.o binary_tree.o parameters.o reaction.o observable.o nucleus.o cortex.o simulation.o

# Compiler and executable name
EXE = gillespie.X$(OPT64)
CFLAGS = $(GFLAG)
GFLAG = #-g3
LDFLAGS=
COMP = g++ -O3 $(CFLAGS) $(LDFLAGS)

# Archive file 'made' remembering the last compilation
MADE = ${strip $(shell cat made)}

ifeq ($(MADE),cluster)
	OPT64 = 64
endif

normal : OPT64 =
normal : REMEMBER = normal
normal : $(EXE)
cluster: CFLAGS = $(GFLAG) -m64
cluster: OPT64 = 64
cluster: REMEMBER = cluster
cluster : $(EXE)
$(EXE) : $(OF) gillespie.o analyse_trajectory.o temp2init.o
	$(COMP) -o $(EXE) gillespie.o $(OF) -lm
	@echo "  Made $(EXE)"
	$(COMP) -o temp2init.X$(OPT64) temp2init.o $(OF) -lm
	@echo "  Made temp2init.X$(OPT64)"
	$(COMP) -o analyse_trajectory.X$(OPT64) analyse_trajectory.o histo1D.o histo2D.o random_variable.o -lm
	@echo "  Made analyse_trajectory.X$(OPT64)"
	$(shell echo $(REMEMBER) > made)
random_variable.o : $(TOOLS)random_variable.cpp $(TOOLS)random_variable.hpp makefile
	$(COMP) -c $(TOOLS)random_variable.cpp 
histo1D.o : $(TOOLS)histo1D.cpp $(TOOLS)histo1D.hpp makefile
	$(COMP) -c $(TOOLS)histo1D.cpp 
histo2D.o : $(TOOLS)histo2D.cpp $(TOOLS)histo2D.hpp $(TOOLS)array2D.hpp makefile
	$(COMP) -c $(TOOLS)histo2D.cpp 
binary_tree.o : $(TOOLS)binary_tree.cpp $(TOOLS)binary_tree.hpp $(TOOLS)tree_node.hpp makefile
	$(COMP) -c $(TOOLS)binary_tree.cpp 
parameters.o : parameters.cpp parameters.hpp observable.hpp observable.cpp makefile
	$(COMP) -c parameters.cpp 
reaction.o : reaction.hpp
	$(COMP) -c reaction.cpp 
observable.o : observable.hpp
	$(COMP) -c observable.cpp 
nucleus.o : nucleus.cpp nucleus.hpp reaction.hpp observable.hpp parameters.hpp $(TOOLS)array2D.hpp $(TOOLS)mersenne_twister.hpp makefile
	$(COMP) -c nucleus.cpp 
cortex.o : cortex.cpp cortex.hpp nucleus.hpp parameters.hpp $(TOOLS)binary_tree.hpp $(TOOLS)mersenne_twister.hpp makefile
	$(COMP) -c cortex.cpp 
simulation.o : simulation.cpp simulation.hpp cortex.hpp parameters.hpp $(TOOLS)random_variable.hpp $(TOOLS)histo1D.hpp $(TOOLS)mersenne_twister.hpp makefile
	$(COMP) -c simulation.cpp
gillespie.o : gillespie.cpp parameters.hpp simulation.hpp makefile
	$(COMP) -c gillespie.cpp
analyse_trajectory.o : analyse_trajectory.cpp analyse_trajectory.hpp histo1D.o histo2D.o random_variable.o $(TOOLS)array2D.hpp $(TOOLS)mersenne_twister.hpp $(TOOLS)random_variable.hpp makefile
	$(COMP) -c analyse_trajectory.cpp
analyse_trajectory : analyse_trajectory.o
	$(COMP) -o analyse_trajectory.X analyse_trajectory.o histo1D.o histo2D.o random_variable.o -lm
temp2init.o : temp2init.cpp makefile
	$(COMP) -c temp2init.cpp
temp2init : temp2init.o
	$(COMP) -o temp2init.X temp2init.o run_parameter.o -lm
calc_corr : calc_corr.cpp
	$(COMP) -o calc_corr.X calc_corr.cpp
clean :
	rm $(OF)
	$(shell echo "clean" > made)
clear :
	rm $(OF)
	$(shell echo "clear" > made)
virgin :
	rm -v *.X *.X64 *.o
	$(shell echo "virgin" > made)
nothing: ;
