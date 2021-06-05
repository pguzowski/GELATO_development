file(GLOB GDCGCore_sources "src/*.cpp")
add_library(GDCGCore INTERFACE )
target_link_libraries(GDCGCore INTERFACE GDCGCore)
