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
CMAKE_SOURCE_DIR = /Users/neerajsarna/sciebo/DG_moments_deal/Tests

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/neerajsarna/sciebo/DG_moments_deal/Tests/build

# Include any dependencies generated for this target.
include CMakeFiles/test.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test.out.dir/flags.make

CMakeFiles/test.out.dir/main.cc.o: CMakeFiles/test.out.dir/flags.make
CMakeFiles/test.out.dir/main.cc.o: ../main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/neerajsarna/sciebo/DG_moments_deal/Tests/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test.out.dir/main.cc.o"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test.out.dir/main.cc.o -c /Users/neerajsarna/sciebo/DG_moments_deal/Tests/main.cc

CMakeFiles/test.out.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.out.dir/main.cc.i"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/neerajsarna/sciebo/DG_moments_deal/Tests/main.cc > CMakeFiles/test.out.dir/main.cc.i

CMakeFiles/test.out.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.out.dir/main.cc.s"
	/Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/neerajsarna/sciebo/DG_moments_deal/Tests/main.cc -o CMakeFiles/test.out.dir/main.cc.s

CMakeFiles/test.out.dir/main.cc.o.requires:

.PHONY : CMakeFiles/test.out.dir/main.cc.o.requires

CMakeFiles/test.out.dir/main.cc.o.provides: CMakeFiles/test.out.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/test.out.dir/build.make CMakeFiles/test.out.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/test.out.dir/main.cc.o.provides

CMakeFiles/test.out.dir/main.cc.o.provides.build: CMakeFiles/test.out.dir/main.cc.o


# Object files for target test.out
test_out_OBJECTS = \
"CMakeFiles/test.out.dir/main.cc.o"

# External object files for target test.out
test_out_EXTERNAL_OBJECTS =

test.out: CMakeFiles/test.out.dir/main.cc.o
test.out: CMakeFiles/test.out.dir/build.make
test.out: /Applications/deal.II.app/Contents/Resources/lib/libdeal_II.8.4.1.dylib
test.out: /usr/local/lib/libgtest.a
test.out: /usr/lib/libbz2.dylib
test.out: /usr/lib/libz.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libparmetis.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libmetis.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtrilinoscouplings.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpiro.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librol.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_muelu.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_ifpack2.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_amesos2.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_tpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos_sacado.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstokhos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librythmos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu-adapters.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu-interface.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libmuelu.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocathyra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocaepetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liblocalapack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libloca.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnoxepetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnoxlapack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libnox.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libintrepid.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteko.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosbelos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosaztecoo.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosamesos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosml.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libstratimikosifpack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack2-adapters.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack2.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasazitpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libModeLaplace.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasaziepetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libanasazi.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkomplex.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libamesos2.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelostpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelosepetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libbelos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libml.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libifpack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libzoltan2.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpamgen_extras.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libpamgen.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libamesos.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libgaleri-xpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libgaleri-epetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libaztecoo.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libisorropia.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/liboptipack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libxpetra-sup.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libxpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyratpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyraepetraext.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyraepetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libthyracore.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libepetraext.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraext.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetrainout.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkostsqr.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetrakernels.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassiclinalg.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassicnodeapi.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtpetraclassic.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libtriutils.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libglobipack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libshards.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libzoltan.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libepetra.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libsacado.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/librtop.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoskokkoscomm.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoskokkoscompat.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosremainder.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosnumerics.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoscomm.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchosparameterlist.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libteuchoscore.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkosalgorithms.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkoscontainers.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5/lib/libkokkoscore.dylib
test.out: /usr/lib/liblapack.dylib
test.out: /usr/lib/libblas.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/lib/libmpi_cxx.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/arpack-ng-d66b8b4/lib/libparpack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/arpack-ng-d66b8b4/lib/libarpack.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libhdf5_hl.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libhdf5.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBO.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBool.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKBRep.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKernel.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKFeat.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKFillet.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKG2d.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKG3d.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKGeomAlgo.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKGeomBase.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKHLR.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKIGES.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKMath.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKMesh.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKOffset.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKPrim.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKShHealing.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEP.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEPAttr.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTEPBase.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKSTL.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKTopAlgo.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/oce-c0cdb53/lib/libTKXSBase.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/p4est-8d811a8/lib/libp4est.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/p4est-8d811a8/lib/libsc.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/openmpi-1.10.2/lib/libmpi.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/slepc-2c065dd/lib/libslepc.dylib
test.out: /Applications/deal.II.app/Contents/Resources/opt/petsc-3e25e16/lib/libpetsc.dylib
test.out: CMakeFiles/test.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/neerajsarna/sciebo/DG_moments_deal/Tests/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test.out.dir/build: test.out

.PHONY : CMakeFiles/test.out.dir/build

CMakeFiles/test.out.dir/requires: CMakeFiles/test.out.dir/main.cc.o.requires

.PHONY : CMakeFiles/test.out.dir/requires

CMakeFiles/test.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test.out.dir/clean

CMakeFiles/test.out.dir/depend:
	cd /Users/neerajsarna/sciebo/DG_moments_deal/Tests/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/neerajsarna/sciebo/DG_moments_deal/Tests /Users/neerajsarna/sciebo/DG_moments_deal/Tests /Users/neerajsarna/sciebo/DG_moments_deal/Tests/build /Users/neerajsarna/sciebo/DG_moments_deal/Tests/build /Users/neerajsarna/sciebo/DG_moments_deal/Tests/build/CMakeFiles/test.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test.out.dir/depend

