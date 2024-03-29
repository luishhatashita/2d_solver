file(GLOB GRID_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

add_library(grid ${GRID_FILES})

#target_sources(2DSolver
#    PRIVATE
#    grid_files
#)

target_include_directories(grid PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(grid PUBLIC io malloc)
