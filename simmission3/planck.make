PKG:=simmission3

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE_F+= $(F_MODPATH)$(OD)

LIB_$(PKG):=$(LIBDIR)/libsimmission.a
LIBOBJ:=general_const.o general_error.o \
	general_maths.o general_matrix.o general_rand.o \
	general_time.o general_vector.o levels_output.o \
	planck_analyticaldynamics.o planck_idealdynamics.o \
	planck_io.o planck_l2orbit.o planck_mission.o \
	planck_mission_new.o planck_missionio.o planck_missionio_new.o \
	planck_pointing.o planck_pointingperiod.o planck_satellite.o \
	planck_scanning.o planck_scanning_new.o simmission3_main.o \
	simmission4_main.o solarsystem_l2.o solarsystem_l2orbit.o \
	solarsystem_orbit.o solarsystem_planet.o solarsystem_solsys.o \
	solarsystem_star.o
ALLOBJ:=$(LIBOBJ) simmission3.o  simmission4.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)
F90BIN:=simmission3 simmission4
MODULES:=simmission3_main simmission4_main
MODULES:=$(MODULES:%=$(OD)/%.$(MOD))

ODEP:=$(LIB_Modules) $(LIB_portability)
$(LIB_$(PKG)): $(LIBOBJ)
$(ALLOBJ): $(ODEP) | $(OD)_mkdir

BDEP:=$(LIB_$(PKG)) $(LIB_Modules) $(LIB_portability) $(LIB_libcfitsio)

$(OD)/simmission3_main.o: $(OD)/general_rand.o $(OD)/general_error.o $(OD)/solarsystem_solsys.o $(OD)/planck_mission.o $(OD)/planck_io.o $(OD)/levels_output.o
$(OD)/simmission3.o: $(OD)/simmission3_main.o
$(OD)/simmission4_main.o: $(OD)/general_rand.o $(OD)/general_error.o $(OD)/planck_mission_new.o $(OD)/planck_io.o $(OD)/levels_output.o
$(OD)/simmission4.o: $(OD)/simmission4_main.o

$(OD)/general_maths.o: $(OD)/general_const.o
$(OD)/general_matrix.o: $(OD)/general_error.o
$(OD)/general_rand.o: $(OD)/general_error.o
$(OD)/general_time.o: $(OD)/general_const.o $(OD)/general_error.o
$(OD)/general_vector.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_maths.o

$(OD)/levels_output.o: $(OD)/general_const.o $(OD)/general_time.o $(OD)/general_vector.o

$(OD)/solarsystem_l2.o: $(OD)/general_error.o
$(OD)/solarsystem_l2orbit.o: $(OD)/general_time.o $(OD)/general_vector.o $(OD)/solarsystem_l2.o $(OD)/solarsystem_orbit.o $(OD)/solarsystem_planet.o
$(OD)/solarsystem_orbit.o: $(OD)/general_error.o $(OD)/general_time.o $(OD)/general_vector.o
$(OD)/solarsystem_planet.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_time.o $(OD)/general_vector.o $(OD)/solarsystem_l2.o $(OD)/solarsystem_orbit.o
$(OD)/solarsystem_solsys.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_time.o $(OD)/solarsystem_l2orbit.o $(OD)/solarsystem_orbit.o $(OD)/solarsystem_planet.o $(OD)/solarsystem_star.o
$(OD)/solarsystem_star.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o

$(OD)/planck_analyticaldynamics.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_maths.o $(OD)/planck_idealdynamics.o
$(OD)/planck_idealdynamics.o: $(OD)/general_vector.o
$(OD)/planck_io.o:
$(OD)/planck_l2orbit.o: $(OD)/general_const.o $(OD)/general_time.o $(OD)/solarsystem_planet.o $(OD)/solarsystem_l2orbit.o
$(OD)/planck_mission.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_time.o $(OD)/general_matrix.o $(OD)/solarsystem_l2orbit.o $(OD)/solarsystem_solsys.o $(OD)/solarsystem_orbit.o $(OD)/solarsystem_planet.o $(OD)/solarsystem_star.o $(OD)/levels_output.o $(OD)/planck_analyticaldynamics.o $(OD)/planck_idealdynamics.o $(OD)/planck_io.o $(OD)/planck_l2orbit.o $(OD)/planck_pointing.o $(OD)/planck_pointingperiod.o $(OD)/planck_satellite.o $(OD)/planck_scanning.o
$(OD)/planck_mission_new.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_time.o $(OD)/solarsystem_star.o $(OD)/levels_output.o $(OD)/planck_analyticaldynamics.o $(OD)/planck_idealdynamics.o $(OD)/planck_io.o $(OD)/planck_pointing.o $(OD)/planck_pointingperiod.o $(OD)/planck_satellite.o $(OD)/planck_scanning_new.o
$(OD)/planck_missionio.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_time.o
$(OD)/planck_missionio_new.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_time.o
$(OD)/planck_pointing.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_rand.o
$(OD)/planck_pointingperiod.o: $(OD)/general_const.o $(OD)/general_vector.o $(OD)/solarsystem_star.o $(OD)/planck_analyticaldynamics.o $(OD)/planck_satellite.o
$(OD)/planck_satellite.o: $(OD)/general_vector.o $(OD)/general_matrix.o
$(OD)/planck_scanning.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_time.o $(OD)/solarsystem_planet.o $(OD)/solarsystem_solsys.o $(OD)/solarsystem_orbit.o $(OD)/planck_missionio.o
$(OD)/planck_scanning_new.o: $(OD)/general_const.o $(OD)/general_error.o $(OD)/general_vector.o $(OD)/general_time.o $(OD)/planck_missionio_new.o

F90BIN:=$(F90BIN:%=$(BINDIR)/%)
$(F90BIN): $(BINDIR)/%: $(OD)/%.o $(BDEP)

all_lib+=$(LIB_$(PKG))
all_f90bin+=$(F90BIN)
all_mod+=$(MODULES)

levels_bin+=$(F90BIN)
