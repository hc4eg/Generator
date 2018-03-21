cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(generate)


IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


if( DEFINED ${BH_cross_sections_dir} )
	message( "BH_cross_sections_dir is defined" )
else( DEFINED ${BH_cross_sections_dir} )
	set( BH_cross_sections_dir "${PROJECT_SOURCE_DIR}/../BH_cross_sections" )
endif( DEFINED ${BH_cross_sections_dir} )

# Include directories need for this compile
include_directories("${BH_cross_sections_dir}")

message("BH_cross_sections_dir = ${BH_cross_sections_dir}" )

# get the extra libraries needed
find_library(BH_cross_sections_lib BH_cross_sections "${BH_cross_sections_dir}" )

# Add the executable
add_executable(generate generate.cpp) 
target_link_libraries(generate ${BH_cross_sections_lib})
 
