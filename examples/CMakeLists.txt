add_executable(example1 example1.cpp)
target_compile_features(example1 PUBLIC cxx_std_14)
target_link_libraries(example1 PRIVATE dkgen)

add_executable(example2 example2.cpp)
target_compile_features(example2 PUBLIC cxx_std_14)
target_link_libraries(example2 PRIVATE dkgen dkgen-physics)

find_package(GSL)
find_package(ROOT)
add_executable(example3 example3.cpp)
target_compile_features(example3 PUBLIC cxx_std_14)
target_link_libraries(example3 PRIVATE dkgen dkgen-physics GSL::gsl ROOT::Tree ROOT::Physics)
