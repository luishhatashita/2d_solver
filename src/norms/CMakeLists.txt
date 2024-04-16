file(GLOB NORMS_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

add_library(norms ${NORMS_FILES})

#target_sources(2DSolver
#    PRIVATE
#    grid_files
#)

target_include_directories(norms PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#target_link_libraries(io PUBLIC else)
