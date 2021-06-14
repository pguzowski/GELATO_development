cmake_minimum_required(VERSION 3.8)

project(dkgen VERSION 1.0.0 LANGUAGES CXX)

add_compile_options(-Wall -Wextra -pedantic -Werror)
if(CMAKE_BUILD_TYPE=Release)
  add_compile_options(-O3)
elseif(CMAKE_BUILD_TYPE=Debug)
  add_compile_options(-O0 -g -DDEBUG)
endif()


add_subdirectory(core)
add_subdirectory(physics)
add_subdirectory(examples)
