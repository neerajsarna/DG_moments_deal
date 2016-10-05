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


SET(DEAL_II_DIR "/Applications/deal.II.app/Contents/Resources/lib/cmake")
SET(MKL_DIR "/opt/intel//compilers_and_libraries_2016.3.170/mac/mkl")

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
#adding the executables

ADD_EXECUTABLE(test.out main.cc)


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../include_new)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include_tests_systemA)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../periodic_square)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../systems)
INCLUDE_DIRECTORIES(/usr/local/include/eigen)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../ring)

# Locate GTest
FIND_PACKAGE(GTest REQUIRED)
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(test.out ${GTEST_LIBRARIES} libmkl_intel_ilp64.dylib libmkl_core.dylib libmkl_intel_thread.dylib libpthread.dylib libdl.dylib pthread)

#TARGET_LINK_LIBRARIES(test.out ${GTEST_LIBRARIES} pthread)

DEAL_II_SETUP_TARGET(test.out)
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

