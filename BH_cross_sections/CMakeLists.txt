# This CMake just creates the cross_sections library
# After running cmake just type make
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(BH_cross_sections)

#IF(NOT CMAKE_BUILD_TYPE)
#  SET(CMAKE_BUILD_TYPE Debug CACHE STRING
#      "Choose the type of build, options are: None Debug Release
#RelWithDebInfo MinSizeRel."
#      FORCE)
#ENDIF(NOT CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE None)

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O0")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
# Unfortunately compiling with anything other than no
# optimizing causes Numerical Recipies to give run
# time errors!

include_directories(${PROJECT_SOURCE_DIR}/.)


file(GLOB xsecheaders ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB xsecsources ${PROJECT_SOURCE_DIR}/*.cpp)
add_library(BH_cross_sections SHARED ${xsecheaders} ${xsecsources})

#install( TARGETS LightOutput
#	LIBRARY DESTINATION lib )
#install( FILES include/LightOutput.hh DESTINATION include )
