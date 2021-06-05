add_executable(test1 test1.cpp)
target_compile_features(test1 PUBLIC cxx_std_14)
target_link_libraries(test1 PRIVATE dkgen)

add_executable(test2 test2.cpp)
target_compile_features(test2 PUBLIC cxx_std_14)
target_link_libraries(test2 PRIVATE dkgen dkgen-physics)
