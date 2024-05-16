cmake_minimum_required(VERSION 3.0)
project(nagata)
FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
FIND_PACKAGE( pybind11 REQUIRED )

FIND_PACKAGE(PythonInterp 3)
FIND_PACKAGE(PythonLibs   3)

IF( PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND AND PYBIND11_FOUND )
  INCLUDE_DIRECTORIES(
    ${PYTHON_INCLUDE_DIRS}
    ${PYBIND11_INCLUDE_DIRS}
  )
  ENDIF()
INCLUDE_DIRECTORIES("/usr/include/python3.7m")

add_library(pynagata SHARED pythonlib.cpp)
 # The library must not have any prefix and should be located in
  # a subfolder that includes the package name. The setup will be
  # more complicated otherwise.
  SET_TARGET_PROPERTIES( pynagata
    PROPERTIES
      PREFIX ""
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pynagata"
  )

add_executable(nagata main.cpp)
install(TARGETS nagata RUNTIME DESTINATION bin)
