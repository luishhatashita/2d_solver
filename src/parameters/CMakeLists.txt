file(GLOB PAR_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

add_library(parameters ${PAR_FILES})

#target_sources(2DSolver
#    PRIVATE
#    grid_files
#)

set_target_properties(parameters PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(parameters PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#target_link_libraries(parameters PUBLIC schemes)
