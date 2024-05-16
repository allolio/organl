cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
#set (STATIC_LINK 1)
if(STATIC_LINK)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
project(organl)
FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
message("OPENMP FOUND")
        # set(PYTHON_INCLUDE_DIR "/usr/include/python3.6")
        # set(PYTHON_INCLUDE_DIRS "/usr/include/python3.6")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#SET(PYTHON_LIBRARIES
#"/usr/lib/python3.10/config-3.6m-x86_64-linux-gnu/libpython3.6m.so")
#SET(PYTHON_INCLUDE_DIRS "/usr/include/python3.10")
#INCLUDE_DIRECTORIES("/usr/include/python3.8")
#FIND_PACKAGE( pybind11 REQUIRED )

#FIND_PACKAGE(PythonInterp 3)
#FIND_PACKAGE(PythonLibs   3)

#IF( PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND AND PYBIND11_FOUND )
        # INCLUDE_DIRECTORIES(
          #    ${PYTHON_INCLUDE_DIRS}
          # ${PYBIND11_INCLUDE_DIRS}
          # )#
          #  ENDIF()
#INCLUDE_DIRECTORIES("/usr/include/python3.10")

#add_library(pynagata SHARED pythonlib.cpp)
 # The library must not have any prefix and should be located in
  # a subfolder that includes the package name. The setup will be
  # more complicated otherwise.
  #   SET_TARGET_PROPERTIES( pynagata
  #   PROPERTIES

  #    PREFIX ""
  #  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pynagata"
        #  )
add_executable(organl main.cpp)
set_target_properties(organl PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(organl PROPERTIES LINK_SEARCH_END_STATIC 1)
#add_executable(nag_mix main_mix.cpp)
add_executable(angle angle.cpp)
add_executable(helfins helfins.cpp)
install(TARGETS organl RUNTIME DESTINATION bin)
