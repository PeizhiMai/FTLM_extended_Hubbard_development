cmake_minimum_required(VERSION 3.16)
project(ftlm_extended_hubbard LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FTLM_LIBOMP_ROOT "" CACHE PATH
  "Optional root for libomp when building with Apple Clang (for example /opt/homebrew/opt/libomp)")

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND FTLM_LIBOMP_ROOT)
  set(OpenMP_C_FLAGS "-Xclang -fopenmp")
  set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "${FTLM_LIBOMP_ROOT}/lib/libomp.dylib")
  set(OpenMP_C_INCLUDE_DIR "${FTLM_LIBOMP_ROOT}/include")
  set(OpenMP_CXX_INCLUDE_DIR "${FTLM_LIBOMP_ROOT}/include")
endif()

add_executable(ftlm_n_vs_mu
  src/main.cpp
  src/thermo_checkpoint.cpp
)

find_package(LAPACK REQUIRED)
find_package(OpenMP)

target_link_libraries(ftlm_n_vs_mu PRIVATE LAPACK::LAPACK)
if(OpenMP_CXX_FOUND)
  target_link_libraries(ftlm_n_vs_mu PRIVATE OpenMP::OpenMP_CXX)
  message(STATUS "FTLM OpenMP enabled")
else()
  message(STATUS "FTLM OpenMP not found; building serial FTLM executable")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  target_compile_options(ftlm_n_vs_mu PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

add_executable(ftlm_ed_n_vs_mu
  src/ed_main.cpp
  src/thermo_checkpoint.cpp
)

target_link_libraries(ftlm_ed_n_vs_mu PRIVATE LAPACK::LAPACK)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  target_compile_options(ftlm_ed_n_vs_mu PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

add_executable(ftlm_ed_conductivity
  src/cond_ed_main.cpp
  src/thermo_checkpoint.cpp
)

target_link_libraries(ftlm_ed_conductivity PRIVATE LAPACK::LAPACK)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  target_compile_options(ftlm_ed_conductivity PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

add_executable(ftlm_conductivity
  src/cond_ftlm_main.cpp
  src/thermo_checkpoint.cpp
)

target_link_libraries(ftlm_conductivity PRIVATE LAPACK::LAPACK)
if(OpenMP_CXX_FOUND)
  target_link_libraries(ftlm_conductivity PRIVATE OpenMP::OpenMP_CXX)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  target_compile_options(ftlm_conductivity PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

add_executable(ftlm_reduce_checkpoint
  src/reduce_main.cpp
  src/thermo_checkpoint.cpp
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  target_compile_options(ftlm_reduce_checkpoint PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(
    NAME incremental_checkpoint
    COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_incremental_checkpoint.py"
      --ftlm "$<TARGET_FILE:ftlm_n_vs_mu>"
      --reducer "$<TARGET_FILE:ftlm_reduce_checkpoint>"
  )
  set_tests_properties(incremental_checkpoint PROPERTIES TIMEOUT 120)
endif()
