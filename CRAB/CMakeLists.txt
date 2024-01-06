#----------------------------------------------------------------------------
# Setup the project
#cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_minimum_required(VERSION 3.16...3.21)
project(CRAB)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

option(With_MultiThread "WITH_MultiThread" OFF)
option(WITH_SYS_TETGEN "WITH_TETGEN" OFF)
option(With_Opticks OFF)


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
#----------------------------------------------------------------------------
# Setup include directory for this project
#

include_directories ("./include/CRAB/")
include(${Geant4_USE_FILE})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT QUIET REQUIRED COMPONENTS RIO Net)
if(NOT ROOT_FOUND)
    message(STATUS "ROOT package not found.")
    return()
endif()


## Garfield
if(NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
endif()

if(NOT Garfield_FOUND)
    message(STATUS "Garfield package not found.")
    return()
endif()

## NEST
find_package(NEST REQUIRED)
if(NOT NEST_FOUND)
    message(STATUS "NEST package not found.")
    return()
endif()

if(With_Opticks)

    message(STATUS "WITH_OPTICKS is set ON")
    list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
    #list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")
    set(OPTICKS_PREFIX "$ENV{OPTICKS_PREFIX}")
    ### we need to Make sure this OPTICKS_PREFIX matches with the one on OpticksBuildOptions
    message(STATUS "OPTICKS_PREFIX : $ENV{OPTICKS_PREFIX} ")
    include(OpticksBuildOptions)
    message(STATUS "---- OPTICKS ----")
    message(STATUS "CMAKE_MODULE_PATH --> ${CMAKE_MODULE_PATH}")
    message(STATUS "OPTICKS HOME --> $ENV{OPTICKS_HOME}")
    message(STATUS "---- OPTICKS END ----")
    find_dependency(PLog) # For Some Reason it was not finding PLog , doing this compiles
    ## These are needed for opticks
    find_package(G4CX REQUIRED CONFIG)
    find_package(U4 REQUIRED CONFIG)
    find_package(ExtG4 REQUIRED CONFIG)
    find_package(CSG_GGeo REQUIRED CONFIG)


endif ()

# Setup Geant4 include directories and compile definitions
include_directories($ENV{PROJECT_SOURCE_DIR}/include)



## Include tetgen
if(${WITH_SYS_TETGEN} MATCHES "ON")
    find_package(tetgen)
endif()


IF(APPLE)
    # Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
    LINK_DIRECTORIES(/usr/local/lib)
ENDIF()

#--------------------------------------------------------------------------
#Create a shared library with the generated dictionary

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/*/*.hh)

# Copy header files to the 'include' binary folder
file(COPY ${headers} DESTINATION ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME})



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${NEST_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${Geant4_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and GARFIELD++ libraries
#

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g")
add_executable(CRAB CRAB.cc ${sources} ${headers})
message(STATUS "GEANT4 libraries libs: ${Geant4_LIBRARIES}")
message(STATUS "NEST Includes: ${NEST_INCLUDE_DIRS}")

## makes finding fortran little easier
#find_library(GFORTRAN_LIBRARY NAMES libgfortran)

target_link_libraries(CRAB -lgfortran Garfield::Garfield ROOT::RIO ROOT::Net NEST::NESTCore NEST::NESTG4 ${Geant4_LIBRARIES} )

if(${WITH_SYS_TETGEN} MATCHES "ON")
    target_link_libraries(CRAB tet)
endif()
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#####################################
if (With_Opticks)
    target_link_libraries(CRAB Opticks::G4CX Opticks::U4 Opticks::ExtG4 Opticks::CSG_GGeo)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
endif()
#Configure some of the variables
configure_file(config.h.in ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/config.h)

#####################################


#target_link_libraries(CRAB -lGarfield -lgfortran ${ROOT_LIBRARIES} libNESTCore.dylib  libNESTG4.dylib ${Geant4_LIBRARIES})
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(CRAB_SCRIPTS
  #Xenon.gas
  #vis.mac
  #run1.mac
  #convertDegradFile.py
)

foreach(_script ${CRAB_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CRAB DESTINATION bin)
