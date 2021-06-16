add_library(Physics INTERFACE )
target_include_directories(Physics INTERFACE 
  #$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  #$<INSTALL_INTERFACE:include>
  include
)
target_link_libraries(Physics INTERFACE dkgen)
target_compile_features(Physics INTERFACE cxx_std_14)

install(TARGETS Physics
  INCLUDES DESTINATION include
  )
install(DIRECTORY ${PROJECT_SOURCE_DIR}/physics/include/dkgen DESTINATION include)
