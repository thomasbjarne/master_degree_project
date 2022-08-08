#finally it works...
FC = gfortran -fcheck=all -g -fbacktrace

main 	= $(wildcard *.f90)
program = $(patsubst %.f90, %, $(main))
mods 	= $(wildcard mod/mod*.f90)
objs 	= $(patsubst mod/%.f90, %.o, $(mods))
temp 	= $(patsubst %.o, obj/%.o, $(objs))

.phony=clean debug all clean_data

all : $(program)

$(program) : $(objs)
	@mv $(objs) obj
	@$(FC) $(temp) $(program).o -o $@
	@mv $(program).o obj

$(objs) : $(mods) $(main)
	@$(FC) $^ -Jmod -c

clean :
	@rm -f obj/*.mod obj/*.o mod/*.mod mod/*.o *.mod *.o $(program)

clean_data:
	@rm -f data/*

debug :
	@echo 'MAIN =' $(main)
	@echo 'PROGRAM ='$(program)
	@echo 'MOD = '$(mods)
	@echo 'OBJS = '$(objs)
	@echo 'TEMP = '$(temp)