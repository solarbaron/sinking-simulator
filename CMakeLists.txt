cmake_minimum_required(VERSION 3.24)
project(shipsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_library(shipsim_engine
  engine/core/geometry.cpp
  engine/sim/ship.cpp
)
target_include_directories(shipsim_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(shipsim_engine PRIVATE -Wall -Wextra -Wpedantic)

add_executable(shipsim
  game/prototype/main.cpp
  game/prototype/ferry.cpp
)
target_link_libraries(shipsim PRIVATE shipsim_engine)
target_compile_options(shipsim PRIVATE -Wall -Wextra)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()
add_executable(shipsim_tests tests/test_core.cpp)
target_link_libraries(shipsim_tests PRIVATE shipsim_engine)
target_compile_options(shipsim_tests PRIVATE -Wall -Wextra)
add_test(NAME core COMMAND shipsim_tests)
