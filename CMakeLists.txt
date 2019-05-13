#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#add_subdirectory(Xenon)
add_subdirectory(ALICE)
add_subdirectory(MWGPC)

#add_custom_target(Interface DEPENDS Xenon ALICE)
add_custom_target(Interface DEPENDS ALICE MWGPC)

