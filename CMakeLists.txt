cmake_minimum_required(VERSION 3.15)

# The kind of compilation : Release or Debug ######################################
set ( CMAKE_BUILD_TYPE "Release" CACHE STRING "Debug or Release" FORCE)
# RelWithDebInfo    Release   Debug

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

project(Shallow
LANGUAGES Fortran
VERSION 1.0
DESCRIPTION "Finite volume solver for the 2D shallow water equations")

include( "CMake.config" )

# Include the headers for the external libraries ###############################
#include_directories(  )

# Include the flags in the compiler options #####################################
if ( Windows )
   message ( "   Windows environment ..." )
   set ( ENVI_FLAG  "-DWINDOWS" )
else ()
   message ( "   Linux/MAX environment..." )
   set ( ENVI_FLAG  "" )
endif()

add_definitions ( ${ENVI_FLAG} ${SIGWATCH_FLAG} )

# The executables are put in the source_dir/EXE directory ##########################
message ( "The executables are to be put in ${CMAKE_CURRENT_SOURCE_DIR}/EXE" )
set ( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/EXE )
set ( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/EXE )

# The compiled modules are put in the source_dir/OBJ directory #####################
set ( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/OBJ )

# Look for specific files in the SRC directory ####################################
file ( GLOB_RECURSE  SRC_file SRC/ module_shallow.f90 )

# Extract the base of the paths of the previous files #############################
get_filename_component (SRC_dir ${SRC_file} PATH )

# Construction of the executables ###########################################
add_library ( libmodules
				      ${SRC_dir}/module_shallow.f90
)

#
add_executable ( shallow 
				      ${SRC_dir}/booklib.f90
				      ${SRC_dir}/get_normal_to_cell.f90
				      ${SRC_dir}/deltaT.f90
				      ${SRC_dir}/gmsh_operations.f90
				      ${SRC_dir}/find_elem.f90
				      ${SRC_dir}/flux.f90
				      ${SRC_dir}/gnufor.f90
				      ${SRC_dir}/Main.f90
				      ${SRC_dir}/runge_kutta.f90
				      ${SRC_dir}/handling_parameters.f90
			)
add_dependencies( shallow libmodules)
set_target_properties ( shallow PROPERTIES LINKER_LANGUAGE Fortran)

#
add_executable ( build_initial_solution 
				      ${SRC_dir}/Build_initial_condition.f90
				      ${SRC_dir}/gmsh_operations.f90
				      ${SRC_dir}/find_elem.f90
#				      ${SRC_dir}/get_normal_to_cell.f90
				      ${SRC_dir}/handling_parameters.f90
			)
add_dependencies( build_initial_solution libmodules)
set_target_properties ( build_initial_solution PROPERTIES LINKER_LANGUAGE Fortran)

#
add_executable ( get_relative_error_theory 
				      ${SRC_dir}/Get_relative_error_theory.f90
				      ${SRC_dir}/gmsh_operations.f90
				      ${SRC_dir}/find_elem.f90
#				      ${SRC_dir}/get_normal_to_cell.f90
				      ${SRC_dir}/handling_parameters.f90
			)
add_dependencies( get_relative_error_theory libmodules)
set_target_properties ( get_relative_error_theory PROPERTIES LINKER_LANGUAGE Fortran)
            
# Linking with all the external libraries ##########################################
target_link_libraries( shallow	                 libmodules ${BLAS_LIB} ${SIGWATCH_LIBRARIES} )
target_link_libraries( build_initial_solution	 libmodules ${BLAS_LIB} )
target_link_libraries( get_relative_error_theory libmodules ${BLAS_LIB} )