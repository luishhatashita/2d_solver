file(GLOB MALLOC_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

add_library(malloc ${MALLOC_FILES})

#target_sources(2DSolver
#    PRIVATE
#    grid_files
#)

target_include_directories(malloc PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#target_link_libraries(malloc PUBLIC else)
