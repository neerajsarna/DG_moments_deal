# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/neerajsarna/sciebo/DG_moments_deal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/neerajsarna/sciebo/DG_moments_deal/build

# Include any dependencies generated for this target.
include CMakeFiles/dg_moments.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dg_moments.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dg_moments.out.dir/flags.make

CMakeFiles/dg_moments.out.dir/main.cc.o: CMakeFiles/dg_moments.out.dir/flags.make
CMakeFiles/dg_moments.out.dir/main.cc.o: ../main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/neerajsarna/sciebo/DG_moments_deal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dg_moments.out.dir/main.cc.o"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dg_moments.out.dir/main.cc.o -c /Users/neerajsarna/sciebo/DG_moments_deal/main.cc

CMakeFiles/dg_moments.out.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dg_moments.out.dir/main.cc.i"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/neerajsarna/sciebo/DG_moments_deal/main.cc > CMakeFiles/dg_moments.out.dir/main.cc.i

CMakeFiles/dg_moments.out.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dg_moments.out.dir/main.cc.s"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/neerajsarna/sciebo/DG_moments_deal/main.cc -o CMakeFiles/dg_moments.out.dir/main.cc.s

CMakeFiles/dg_moments.out.dir/main.cc.o.requires:

.PHONY : CMakeFiles/dg_moments.out.dir/main.cc.o.requires

CMakeFiles/dg_moments.out.dir/main.cc.o.provides: CMakeFiles/dg_moments.out.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/dg_moments.out.dir/build.make CMakeFiles/dg_moments.out.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/dg_moments.out.dir/main.cc.o.provides

CMakeFiles/dg_moments.out.dir/main.cc.o.provides.build: CMakeFiles/dg_moments.out.dir/main.cc.o


# Object files for target dg_moments.out
dg_moments_out_OBJECTS = \
"CMakeFiles/dg_moments.out.dir/main.cc.o"

# External object files for target dg_moments.out
dg_moments_out_EXTERNAL_OBJECTS =

dg_moments.out: CMakeFiles/dg_moments.out.dir/main.cc.o
dg_moments.out: CMakeFiles/dg_moments.out.dir/build.make
dg_moments.out: /Applications/deal.II.app/Contents/Resources/lib/libdeal_II.g.8.4.1.dylib
dg_moments.out: /usr/lib/libbz2.dylib
dg_moments.out: /usr/lib/libz.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libparmetis.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libmetis.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtrilinoscouplings.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpiro.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librol.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_muelu.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_ifpack2.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_amesos2.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_tpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_sacado.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librythmos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu-adapters.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu-interface.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocathyra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocaepetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocalapack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libloca.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnoxepetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnoxlapack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnox.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libintrepid.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteko.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosbelos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosaztecoo.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosamesos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosml.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosifpack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack2-adapters.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack2.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasazitpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libModeLaplace.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasaziepetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasazi.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkomplex.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libamesos2.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelostpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelosepetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libml.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libzoltan2.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpamgen_extras.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpamgen.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libamesos.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libgaleri-xpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libgaleri-epetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libaztecoo.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libisorropia.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liboptipack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libxpetra-sup.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libxpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyratpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyraepetraext.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyraepetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyracore.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libepetraext.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraext.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetrainout.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkostsqr.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetrakernels.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassiclinalg.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassicnodeapi.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassic.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtriutils.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libglobipack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libshards.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libzoltan.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libepetra.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libsacado.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librtop.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoskokkoscomm.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoskokkoscompat.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosremainder.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosnumerics.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoscomm.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosparameterlist.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoscore.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkosalgorithms.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkoscontainers.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkoscore.dylib
dg_moments.out: /usr/lib/liblapack.dylib
dg_moments.out: /usr/lib/libblas.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/lib/libmpi_cxx.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/arpack-ng-d66b8b4/lib/libparpack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/arpack-ng-d66b8b4/lib/libarpack.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libhdf5_hl.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libhdf5.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBO.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBool.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBRep.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKernel.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKFeat.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKFillet.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKG2d.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKG3d.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKGeomAlgo.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKGeomBase.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKHLR.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKIGES.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKMath.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKMesh.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKOffset.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKPrim.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKShHealing.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEP.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEPAttr.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEPBase.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTL.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKTopAlgo.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKXSBase.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/p4est-8d811a8/lib/libp4est.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/p4est-8d811a8/lib/libsc.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/lib/libmpi.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/slepc-2c065dd/lib/libslepc.dylib
dg_moments.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libpetsc.dylib
dg_moments.out: CMakeFiles/dg_moments.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/neerajsarna/sciebo/DG_moments_deal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dg_moments.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dg_moments.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dg_moments.out.dir/build: dg_moments.out

.PHONY : CMakeFiles/dg_moments.out.dir/build

CMakeFiles/dg_moments.out.dir/requires: CMakeFiles/dg_moments.out.dir/main.cc.o.requires

.PHONY : CMakeFiles/dg_moments.out.dir/requires

CMakeFiles/dg_moments.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dg_moments.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dg_moments.out.dir/clean

CMakeFiles/dg_moments.out.dir/depend:
	cd /Users/neerajsarna/sciebo/DG_moments_deal/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/neerajsarna/sciebo/DG_moments_deal /Users/neerajsarna/sciebo/DG_moments_deal /Users/neerajsarna/sciebo/DG_moments_deal/build /Users/neerajsarna/sciebo/DG_moments_deal/build /Users/neerajsarna/sciebo/DG_moments_deal/build/CMakeFiles/dg_moments.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dg_moments.out.dir/depend

