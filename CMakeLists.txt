cmake_minimum_required(VERSION 3.8)

project(GenericDecayChainGenerator VERSION 1.0.0 LANGUAGES CXX)

find_package(CLHEP REQUIRED)


# test for <experimental/propagate_const> (C++ TS2 standard, not available in all compilers)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx "#include <experimental/propagate_const>\nstd::experimental::propagate_const<int*> x;\nint main(int argc, char** argv) { return 0; }\n")
try_compile(HAS_PROPAGATE_CONST ${CMAKE_CURRENT_BINARY_DIR}/test_compile ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx CXX_STANDARD 14)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx)
file(REMOVE_RECURSE  ${CMAKE_CURRENT_BINARY_DIR}/test_compile)

if(NOT ${HAS_PROPAGATE_CONST})
  message(FATAL_ERROR "Compiler doesn't have std::experimental::propagate_const feature")
endif()

# end test for <experimental/propagate_const>

add_subdirectory(core)
add_subdirectory(test)
