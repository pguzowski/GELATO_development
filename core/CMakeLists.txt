set(PHYSICS_BACKEND "CLHEP" CACHE STRING "Choice of physics backend library, can be CLHEP or ROOT for now. CLHEP is faster than ROOT in benchmarks")

set(EXPOSE_PHYSICS_VECTORS OFF CACHE BOOL "Directly use CLHEP vectors, rotations, etc, instead of hiding these away")

if(EXPOSE_PHYSICS_VECTORS)
  set(PHYSICS_BACKEND "CLHEP")
endif()

if(PHYSICS_BACKEND STREQUAL "CLHEP")
  find_package(CLHEP REQUIRED)
elseif(PHYSICS_BACKEND STREQUAL "ROOT")
  find_package(ROOT REQUIRED)
else()
  message(FATAL_ERROR "PHYSICS_BACKEND variable must be either ROOT or CLHEP")
endif()


if(NOT EXPOSE_PHYSICS_VECTORS)

  #[=[

    # test for <experimental/propagate_const> (C++ TS2 standard, not available in all compilers)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx "#include <experimental/propagate_const>\nstd::experimental::propagate_const<int*> x;\nint main(int argc, char** argv) { return 0; }\n")
    try_compile(RESULT_VAR ${CMAKE_CURRENT_BINARY_DIR}/test_compile ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx CXX_STANDARD 14)
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx)
    file(REMOVE_RECURSE  ${CMAKE_CURRENT_BINARY_DIR}/test_compile)

    if(NOT ${RESULT_VAR})
      message(FATAL_ERROR "Compiler doesn't have std::experimental::propagate_const feature")
    endif()

    # end test for <experimental/propagate_const>

  ]=]#

  # test for physics vector sizes

  function(CHECK_FOR_SIZES)
    if(ARGC GREATER 3)
      #message(STATUS ARGC " = " ${ARGC} " , " ${ARGV0} " , " ${ARGV1} " , " ${ARGV2} " , " ${ARGV3})
      set(INCFILE ${ARGV0})
      set(INCDIRS ${ARGV1})
      set(CLASSNAME ${ARGV2})
      set(OUTVAR ${ARGV3})
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx "#include <${INCFILE}>\nclass derived: public ${CLASSNAME} { };\ntemplate<size_t N> char sizeof_() { return N+256; }\nauto _ = sizeof_<sizeof(derived)>;")
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx "#include \"stddef.h\"\ntemplate<size_t N> char sizeof_() { return N+256; }\nconstexpr size_t M=314159265;\nauto _ = sizeof_<M>;")
    endif()
    try_compile(RESULT_VAR ${CMAKE_CURRENT_BINARY_DIR}/test_compile ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx
      CXX_STANDARD 14
      OUTPUT_VARIABLE RESULT_COMPILE_WORKED
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${INCDIRS}" 
      )
    string(REGEX MATCH
      "(sizeof_<[0-9]+>|[ \t\r\n]N[ \t\r\n]*=[ \t\r\n]*[0-9]+)"
      REGEX_ATTEMPT_RESULT ${RESULT_COMPILE_WORKED})
    string(REGEX MATCH
      "[0-9]+"
      REGEX_ATTEMPT_RESULT_2 ${REGEX_ATTEMPT_RESULT})
    if(ARGC GREATER 3)
      message(STATUS ${OUTVAR} " = '" ${REGEX_ATTEMPT_RESULT_2} "'")
      if(NOT REGEX_ATTEMPT_RESULT_2 STREQUAL "")
        set(${OUTVAR} ${REGEX_ATTEMPT_RESULT_2} PARENT_SCOPE)
      endif()
    else()
      if(NOT REGEX_ATTEMPT_RESULT_2 STREQUAL "314159265")
        message(FATAL_ERROR "Cannot determine sizes of objects at compile time (please contact pawel.guzowski@manchester.ac.uk with more details to attempt a fix)")
      endif()
    endif()
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test_compile.cxx)
    file(REMOVE_RECURSE  ${CMAKE_CURRENT_BINARY_DIR}/test_compile)
  endfunction()

  #first check if regexes will work
  check_for_sizes()
  if(PHYSICS_BACKEND STREQUAL "CLHEP")
    check_for_sizes("CLHEP/Vector/ThreeVector.h" "${CLHEP_INCLUDE_DIRS}" "CLHEP::Hep3Vector" PHYSICS_VECTOR3_IMPL_SIZE)
    check_for_sizes("CLHEP/Vector/LorentzVector.h" "${CLHEP_INCLUDE_DIRS}" "CLHEP::HepLorentzVector" PHYSICS_FOURVECTOR_IMPL_SIZE)
    check_for_sizes("CLHEP/Vector/Rotation.h" "${CLHEP_INCLUDE_DIRS}" "CLHEP::HepRotation" PHYSICS_ROTATION_IMPL_SIZE)
  elseif(PHYSICS_BACKEND STREQUAL "ROOT")
    check_for_sizes("TVector3.h" "${ROOT_INCLUDE_DIRS}" "TVector3" PHYSICS_VECTOR3_IMPL_SIZE)
    check_for_sizes("TLorentzVector.h" "${ROOT_INCLUDE_DIRS}" "TLorentzVector" PHYSICS_FOURVECTOR_IMPL_SIZE)
    check_for_sizes("TRotation.h" "${ROOT_INCLUDE_DIRS}" "TRotation" PHYSICS_ROTATION_IMPL_SIZE)
  endif()

  configure_file(
    "include/GELATO/core/physics_vector_sizes.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/GELATO/core/physics_vector_sizes.hpp"
    )
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/include/GELATO/core/physics_vector_sizes.hpp"
    DESTINATION include/GELATO/core/)
endif()


# end test for physics vector sizes



add_library(Core SHARED
  "src/config.cpp"
  "src/decaying_particle_info.cpp"
  "src/driver.cpp"
  "src/geometry.cpp"
  "src/ordered_list_of_vectors.cpp"
  "src/particle.cpp"
  "src/particle_history.cpp"
  "src/vectors.cpp"
  )
add_library(GELATO::Core ALIAS Core)
target_compile_features(Core PUBLIC cxx_std_14)
if(PHYSICS_BACKEND STREQUAL "CLHEP")
  if(EXPOSE_PHYSICS_VECTORS)
    target_link_libraries(Core PUBLIC CLHEP::Vector)
    target_compile_definitions(Core PUBLIC USING_CLHEP EXPOSE_PHYSICS_VECTORS)
  else()
    target_link_libraries(Core PRIVATE CLHEP::Vector)
    target_compile_definitions(Core PRIVATE USING_CLHEP)
  endif()
elseif(PHYSICS_BACKEND STREQUAL "ROOT")
  target_compile_definitions(Core PRIVATE USING_ROOT)
  target_link_libraries(Core PRIVATE ROOT::Physics)
endif()
target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
if(NOT EXPOSE_PHYSICS_VECTORS)
  target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)
endif()

install(TARGETS Core
  EXPORT GELATOCoreTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )
install(FILES
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/config.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/decaying_particle_info.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/driver.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/geometry.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/ordered_list_of_vectors.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/particle.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/particle_history.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/random_number_generators.hpp"
  "${PROJECT_SOURCE_DIR}/core/include/GELATO/core/vectors.hpp"
  DESTINATION include/GELATO/core/)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  GELATOCoreConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

install(EXPORT GELATOCoreTargets
  FILE GELATOCoreTargets.cmake
  NAMESPACE GELATO::
  DESTINATION lib/cmake/GELATO
  )
