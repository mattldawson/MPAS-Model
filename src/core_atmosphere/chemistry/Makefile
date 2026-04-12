.SUFFIXES: .F .F90 .o

OBJS = mpas_chemistry_utils.o \
       mpas_chemistry_state.o \
       mpas_chemistry_tuvx.o \
       mpas_chemistry_micm.o \
       mpas_chemistry_driver.o

all: $(OBJS)

# Dependencies
mpas_chemistry_state.o:

mpas_chemistry_utils.o:

mpas_chemistry_tuvx.o: mpas_chemistry_state.o

mpas_chemistry_micm.o:

mpas_chemistry_driver.o: mpas_chemistry_micm.o mpas_chemistry_tuvx.o \
                          mpas_chemistry_state.o mpas_chemistry_utils.o

clean:
	$(RM) *.o *.mod *.f90
	$(RM) *.i

.F90.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../../external/esmf_time_f90
endif

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../../external/esmf_time_f90
endif
