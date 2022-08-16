FC = gfortran -fcheck=all -g -fbacktrace

main 	= $(wildcard *.f90)
program = $(patsubst %.f90, %, $(main))
program_obj = $(patsubst %.f90, %.o, $(main))
mods 	= $(wildcard mod/mod*.f90)
objs 	= $(patsubst mod/%.f90, %.o, $(mods))
temp 	= $(patsubst %.o, obj/%.o, $(objs))

.phony=clean debug all clean_data

all : $(program) move

move :
	@mv $(objs) obj
	@mv $(program_obj) obj

$(program) : $(objs)
	@$(FC) $(objs) $@.o -o $@

$(objs) : $(mods) $(main)
	@$(FC) $^ -Jmod -c
	

clean :
	@rm -f obj/*.o mod/*.mod mod/*.o *.mod *.o $(program)
	@rm -f data/data*

debug :
	@echo 'MAIN =' $(main)
	@echo 'PROGRAM ='$(program)
	@echo 'MODS = '$(mods)
	@echo 'OBJS = '$(objs)
	@echo 'TEMP = '$(temp)