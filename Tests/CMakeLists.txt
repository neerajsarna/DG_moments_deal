#setting the (name of the target
SET(TARGET "test")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(G++_COVERAGE_COMPILE_FLAGS "-Wno-unused-local-typedefs")
SET(G++_COVERAGE_COMPILE_FLAGS2 "-pg")
SET(G++_COVERAGE_COMPILE_FLAGS3 "-Wno-pedantic")
SET(G++_COVERAGE_COMPILE_FLAGS4 "-Wno-reorder")
SET(G++_COVERAGE_COMPILE_FLAGS5 "-Wno-missing-field-initializers")
SET(G++_COVERAGE_COMPILE_FLAGS6 "-DMKL_ILP64 -m64")
SET(G++_COVERAGE_COMPILE_FLAGS7 "-Wl,--no-as-needed")
SET(G++_COVERAGE_COMPILE_FLAGS8 "-fmax-errors=2")
SET(G++_COVERAGE_LINK_FLAGS    "-lm")

SET(LAPTOP ON)
SET(TURBINE OFF)
SET(CLUSTER OFF)


IF(LAPTOP)
  MESSAGE("ON LAPTOP")
  SET(DEAL_II_DIR "/Applications/deal.II.app/Contents/Resources/lib/cmake")
  SET(MKL_DIR "/opt/intel//compilers_and_libraries_2016.3.170/mac/mkl")
  SET(EIGEN_DIR "/usr/local/include/eigen")
ENDIF()

IF(TURBINE)
  SET("ON TURBINE")
  SET(DEAL_II_DIR "/home/sarna/Desktop/Libraries/dealii/lib/cmake")
  SET(MKL_DIR "/opt/intel/composer_xe_2013_sp1.0.080/mkl")
  SET(EIGEN_DIR "/home/sarna/Desktop/Libraries/eigenv3/eigen")
ENDIF()


IF(CLUSTER)
  SET("ON CLUSTER")
  SET(DEAL_II_DIR "/home/ns179556/dealii/lib/cmake")
  SET(MKL_DIR "/opt/intel/mkl/2017.0.098/mkl")
  SET(EIGEN_DIR "/home/ns179556/eigen")
  SET(MKL_LIB_DIR "/opt/intel/mkl/10.2.6.038/lib/64")
ENDIF()




FIND_PACKAGE(deal.II 8.3 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()



DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})

ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS2})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS3})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS4})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS5})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS6})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS7})
ADD_DEFINITIONS(${G++_COVERAGE_COMPILE_FLAGS8})
#adding the executables

ADD_EXECUTABLE(test.out main.cc)


SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${G++_COVERAGE_LINK_FLAGS}" )

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../include_new_m_adaptive)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Test_Generic)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../periodic_square)
INCLUDE_DIRECTORIES(${EIGEN_DIR})


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../Projector_Data)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../Symmetrizer_Data)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../systems)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../ring)
INCLUDE_DIRECTORIES(${MKL_DIR}/include)


IF(LAPTOP )
  FIND_PACKAGE(GTest REQUIRED)
ENDIF()

IF(TURBINE)
  SET(GTEST_INCLUDE_DIRS "/home/sarna/Desktop/Libraries/googletest/include")
  SET(GTEST_LIBRARIES "/home/sarna/Desktop/Libraries/googletest/lib/libgtest.a")
ENDIF()


IF(CLUSTER)
  SET(MKL_LIBRARY "/usr/local_rwth/sw/turbomole/tm6.5/libso/em64t-unknown-linux-gnu_mpi/libiomp5.so")
  SET(GTEST_INCLUDE_DIRS "/home/ns179556/google_testing/include")
  SET(GTEST_LIBRARIES "/home/ns179556/google_testing/lib/libgtest.a")
ENDIF()


# Locate GTest
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})

DEAL_II_SETUP_TARGET(test.out)


IF(LAPTOP)
TARGET_LINK_LIBRARIES(test.out ${GTEST_LIBRARIES} libmkl_intel_ilp64.dylib libmkl_core.dylib libmkl_intel_thread.dylib libpthread.dylib libdl.dylib pthread)
ENDIF()

IF(TURBINE)
TARGET_LINK_LIBRARIES(test.out ${GTEST_LIBRARIES} libiomp5.so libpthread.so libmkl_intel_ilp64.so libmkl_core.so libmkl_intel_thread.so libdl.so)
ENDIF()

IF(CLUSTER)
TARGET_LINK_LIBRARIES(test.out ${MKL_LIBRARY} ${GTEST_LIBRARIES} libpthread.so libmkl_intel_ilp64.so libmkl_core.so libmkl_intel_thread.so libdl.so)
ENDIF()


ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )
  
ADD_CUSTOM_TARGET(run COMMAND mycode.out
 COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"  )

