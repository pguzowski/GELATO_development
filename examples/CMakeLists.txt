add_executable(example1 example1.cpp)
target_link_libraries(example1 PRIVATE dkgen)

add_executable(example2 example2.cpp)
target_link_libraries(example2 PRIVATE dkgen Physics)

find_package(GSL)
find_package(ROOT)
add_executable(example3 example3.cpp)

target_link_libraries(example3 PRIVATE dkgen Physics GSL::gsl ROOT::Tree ROOT::Physics)
install(TARGETS example1 example2 example3
        RUNTIME DESTINATION bin
        )
