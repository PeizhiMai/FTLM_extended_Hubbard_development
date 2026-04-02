cmake_minimum_required(VERSION 3.16)
project(ftlm_extended_hubbard LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ftlm_n_vs_mu
  src/main.cpp
)

find_package(LAPACK REQUIRED)

target_link_libraries(ftlm_n_vs_mu PRIVATE LAPACK::LAPACK)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(ftlm_n_vs_mu PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

add_executable(ftlm_ed_n_vs_mu
  src/ed_main.cpp
)

target_link_libraries(ftlm_ed_n_vs_mu PRIVATE LAPACK::LAPACK)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(ftlm_ed_n_vs_mu PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()
