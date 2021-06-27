add_executable(example1 example1.cpp)
target_link_libraries(example1 PRIVATE GELATO::Core)

add_executable(example2 example2.cpp)
target_link_libraries(example2 PRIVATE GELATO::Core GELATO::Physics)

find_package(GSL)
find_package(ROOT)
add_executable(example3 example3.cpp)
target_link_libraries(example3 PRIVATE GELATO::Core GELATO::Physics GSL::gsl ROOT::Tree ROOT::Physics)

#find_package(Boost)
#add_executable(example4 example4.cpp)
#target_link_libraries(example4 PRIVATE Boost::boost )
#target_compile_features(example4 PRIVATE cxx_std_14 )

add_executable(example5 example5.cpp)
target_link_libraries(example5 PRIVATE GELATO::Core GELATO::Physics GSL::gsl  ROOT::Tree ROOT::Physics)


install(TARGETS example1 example2 example3 example5 #example4
        RUNTIME DESTINATION bin
        )
