cmake_minimum_required(VERSION 3.24)
project(shipsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Sanitizer builds. The job system is the reason this exists: a work-stealing
# deque that passes its functional tests can still be riddled with data races,
# and only a sanitizer run will say so.
#   cmake -S . -B build-tsan -DSHIPSIM_SANITIZE=thread
#   cmake -S . -B build-asan -DSHIPSIM_SANITIZE=address
set(SHIPSIM_SANITIZE "" CACHE STRING "Sanitizer to build with: thread, address, or empty")
if(SHIPSIM_SANITIZE STREQUAL "thread")
  add_compile_options(-fsanitize=thread -g -O1)
  add_link_options(-fsanitize=thread)
elseif(SHIPSIM_SANITIZE STREQUAL "address")
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1)
  add_link_options(-fsanitize=address,undefined)
elseif(NOT SHIPSIM_SANITIZE STREQUAL "")
  message(FATAL_ERROR "SHIPSIM_SANITIZE must be 'thread', 'address', or empty")
endif()

add_library(shipsim_engine
  engine/core/arena.cpp
  engine/core/ecs.cpp
  engine/core/geometry.cpp
  engine/core/jobs.cpp
  engine/core/scheduler.cpp
  engine/core/serialise.cpp
  engine/sim/ship.cpp
  engine/sim/fem.cpp
)
target_include_directories(shipsim_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Threads REQUIRED)
target_link_libraries(shipsim_engine PUBLIC Threads::Threads)
target_compile_options(shipsim_engine PRIVATE -Wall -Wextra -Wpedantic)

add_executable(shipsim
  game/prototype/main.cpp
  game/prototype/ferry.cpp
)
target_link_libraries(shipsim PRIVATE shipsim_engine)
target_compile_options(shipsim PRIVATE -Wall -Wextra)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()
add_executable(job_bench tools/job_bench/main.cpp)
target_link_libraries(job_bench PRIVATE shipsim_engine)
target_compile_options(job_bench PRIVATE -Wall -Wextra)

add_executable(shipsim_tests
  tests/main.cpp
  tests/harness.cpp
  tests/test_arena.cpp
  tests/test_core.cpp
  tests/test_ecs.cpp
  tests/test_scheduler.cpp
  tests/test_serialise.cpp
  tests/test_worldio.cpp
  tests/test_jobs.cpp
)
target_include_directories(shipsim_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
target_link_libraries(shipsim_tests PRIVATE shipsim_engine)
target_compile_options(shipsim_tests PRIVATE -Wall -Wextra)
add_test(NAME core COMMAND shipsim_tests)

# --- FEM spike: explicit tetrahedral solver on the GPU -----------------------
# Optional. Without Vulkan or glslc the rest of the project still builds; this is
# a de-risking probe, not a dependency.
find_package(Vulkan QUIET)
find_program(GLSLC glslc)

if(Vulkan_FOUND AND GLSLC)
  set(SHADER_OUT ${CMAKE_BINARY_DIR}/shaders)
  file(MAKE_DIRECTORY ${SHADER_OUT})
  set(SPV_FILES)
  foreach(SHADER tet_forces.comp node_integrate.comp)
    add_custom_command(
      OUTPUT ${SHADER_OUT}/${SHADER}.spv
      COMMAND ${GLSLC} -O
              ${CMAKE_CURRENT_SOURCE_DIR}/engine/gpu/shaders/${SHADER}
              -o ${SHADER_OUT}/${SHADER}.spv
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/engine/gpu/shaders/${SHADER}
      COMMENT "glslc ${SHADER}")
    list(APPEND SPV_FILES ${SHADER_OUT}/${SHADER}.spv)
  endforeach()
  add_custom_target(shipsim_shaders DEPENDS ${SPV_FILES})

  add_executable(fem_spike tools/fem_spike/main.cpp engine/gpu/fem_gpu.cpp)
  add_dependencies(fem_spike shipsim_shaders)
  target_link_libraries(fem_spike PRIVATE shipsim_engine Vulkan::Vulkan)
  target_compile_definitions(fem_spike PRIVATE SHIPSIM_SHADER_DIR="${SHADER_OUT}")
  # The `{VK_STRUCTURE_TYPE_X}` idiom value-initialises the remaining fields, which
  # is exactly what Vulkan wants; the warning is noise here.
  target_compile_options(fem_spike PRIVATE -Wall -Wextra -Wno-missing-field-initializers)
  add_test(NAME fem COMMAND fem_spike)
else()
  message(STATUS "Vulkan or glslc not found - skipping the FEM spike target")
endif()
