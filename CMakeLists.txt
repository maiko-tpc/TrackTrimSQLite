#----------------------------------------------------------------------------
# set up C++11 compiler
#set(CMAKE_CXX_COMPILER "/miho/sharedAPL/gcc-4.9.3-R2/bin/g++")

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(testTrackTrimSQLite)

#----------------------------------------------------------------------------
# Find ROOT package (2018.04.13)
find_package(ROOT REQUIRED COMPONENTS Hist RIO Physics Tree Gpad Graf3d)


#----------------------------------------------------------------------------
# Setup ROOT include directories and compile definitions (2018.04.13)
#
include(${ROOT_USE_FILE})


#----------------------------------------------------------------------------
# Find Garfield++ package (2018.05.22)
#
list(APPEND CMAKE_PREFIX_PATH $ENV{GARFIELD_HOME})
find_package(Garfield)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include
		    ${ROOT_INCLUDE_DIRS}
		    ${GARFIELD_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)


link_directories($ENV{GARFIELD_HOME}/Library)

#----------------------------------------------------------------------------
# Add executables
#




add_executable(makedb makedb.cpp ${sources} ${headers})
target_link_libraries(makedb ${ROOT_LIBRARIES})
target_link_libraries(makedb ${GARFIELD_LIBRARIES})
target_link_libraries(makedb gfortran)
target_link_libraries(makedb sqlite3)
target_compile_options(makedb PRIVATE -std=c++1y)



add_executable(testTrackTrimSQLite testTrackTrimSQLite.cpp ${sources} ${headers})
target_link_libraries(testTrackTrimSQLite ${ROOT_LIBRARIES})
target_link_libraries(testTrackTrimSQLite ${GARFIELD_LIBRARIES})
target_link_libraries(testTrackTrimSQLite gfortran)
target_link_libraries(testTrackTrimSQLite sqlite3)
target_compile_options(testTrackTrimSQLite PRIVATE -std=c++1y)
