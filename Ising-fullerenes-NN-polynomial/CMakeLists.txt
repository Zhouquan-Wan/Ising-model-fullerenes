cmake_minimum_required(VERSION 3.0)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
PROJECT(planar_ising)

set(SRC
${PROJECT_SOURCE_DIR}/main.cpp 
${PROJECT_SOURCE_DIR}/utils.cpp 
)

# set(DEBUGSRC
# ${PROJECT_SOURCE_DIR}/main-debug.cpp 
# ${PROJECT_SOURCE_DIR}/utils.cpp
# )

include_directories(${PROJECT_SOURCE_DIR})

set(BOOST_ROOT "~/libs/boost_1_82_0")
# message(STATUS "${HOME}/libs")

find_package(Boost)
message(STATUS "Boost version: ${Boost_VERSION}")
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(main ${SRC})
  # add_executable(debug ${DEBUGSRC})
endif()