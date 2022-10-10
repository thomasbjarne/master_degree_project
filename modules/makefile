FC = gfortran -fcheck=all -g -fbacktrace

main 	= $(wildcard *.f90)
program = $(patsubst %.f90, %, $(main))
program_obj = $(patsubst %.f90, %.o, $(main))
mods 	= $(wildcard mod/mod*.f90)
objs 	= $(patsubst mod/%.f90, %.o, $(mods))

.phony.=clean debug all move splot pplot

all : $(program) move

move :
	@mv $(objs) obj
	@mv $(program_obj) obj

$(program) : $(objs)
	@$(FC) $(objs) $@.o -o $@

$(objs) : $(mods) $(main)
	@$(FC) $^ -Jmod -c

gplot :
	@gnuplot plt/gridplot.plt
	
pplot :
	@gnuplot plt/pointplot.plt

splot : 
	@gnuplot plt/surfaceplot.plt

clean :
	@rm -f obj/*.o mod/*.o *.o $(program)
	@rm -f data/data*
	@rm -f *.png

debug :
	@echo 'MAIN =' $(main)
	@echo 'PROGRAM ='$(program)
	@echo 'MODS = '$(mods)
	@echo 'OBJS = '$(objs)
	@echo 'TEMP = '$(temp)