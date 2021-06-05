set(PHYSICS_BACKEND "CLHEP" CACHE STRING "Choice of physics backend library, can be CLHEP or ROOT for now. CLHEP is faster than ROOT in benchmarks")

if(PHYSICS_BACKEND STREQUAL "CLHEP")
  find_package(CLHEP REQUIRED)
elseif(PHYSICS_BACKEND STREQUAL "ROOT")
  find_package(ROOT REQUIRED)
else()
  message(FATAL_ERROR "PHYSICS_BACKEND variable must be either ROOT or CLHEP")
endif()


# test for <experimental/propagate_const> (C++ TS2 standard, not available in all compilers)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx "#include <experimental/propagate_const>\nstd::experimental::propagate_const<int*> x;\nint main(int argc, char** argv) { return 0; }\n")
try_compile(HAS_PROPAGATE_CONST ${CMAKE_CURRENT_BINARY_DIR}/test_compile ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx CXX_STANDARD 14)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx)
file(REMOVE_RECURSE  ${CMAKE_CURRENT_BINARY_DIR}/test_compile)

if(NOT ${HAS_PROPAGATE_CONST})
  message(FATAL_ERROR "Compiler doesn't have std::experimental::propagate_const feature")
endif()

# end test for <experimental/propagate_const>



add_library(dkgen SHARED
  "src/config.cpp"
  "src/decaying_particle_info.cpp"
  "src/driver.cpp"
  "src/geometry.cpp"
  "src/ordered_list_of_vectors.cpp"
  "src/particle.cpp"
  "src/particle_history.cpp"
  "src/vectors.cpp"
  )
target_compile_features(dkgen PUBLIC cxx_std_14)
if(PHYSICS_BACKEND STREQUAL "CLHEP")
  target_compile_definitions(dkgen PRIVATE USING_CLHEP)
  target_link_libraries(dkgen PRIVATE CLHEP::Vector)
  target_include_directories(dkgen PUBLIC include PRIVATE ${CLHEP_INCLUDE_DIRS})
elseif(PHYSICS_BACKEND STREQUAL "ROOT")
  target_compile_definitions(dkgen PRIVATE USING_ROOT)
  target_link_libraries(dkgen PRIVATE ROOT::Physics)
  target_include_directories(dkgen PUBLIC include PRIVATE ${ROOT_INCLUDE_DIRS})
endif()
