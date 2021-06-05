add_library(dkgen-physics INTERFACE )
target_include_directories(dkgen-physics INTERFACE include)
target_link_libraries(dkgen-physics INTERFACE dkgen)
