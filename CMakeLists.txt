cmake_minimum_required(VERSION 2.8)
project(DEMPgen)
set (DEMPgen_VERSION_MAJOR 0)
set (DEMPgen_VERSION_MINOR 1)

add_definitions (-DROOT_MATH_LOG)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(include)
include_directories("src/eic_evgen")
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
#include_directories(include/json)

file(GLOB SOURCES "src/*" "src/eic_evgen/eic.cc" "src/eic_evgen/eic_pim.cc" "src/eic_evgen/tssa_sig_Para.cc" "src/eic_evgen/reaction_routine.cc" "src/eic_evgen/legacy_function.cc" "src/eic_evgen/legacy_function.cc" "src/eic_evgen/process_routine/*")
# file(GLOB SOURCES "src/*.cpp")

message("Root Include Dirs: " ${ROOT_INCLUDE_DIRS})
message("Root Libraries: " ${ROOT_LIBRARIES})
message("Mathmore: " ${ROOT_mathmore_FOUND})

add_executable (DEMPgen ${SOURCES})
target_link_libraries(DEMPgen ${ROOT_LIBRARIES})
