file(GLOB SCHEMES_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

add_library(schemes ${SCHEMES_FILES})

#target_sources(2DSolver
#    PRIVATE
#    grid_files
#)

target_include_directories(schemes PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(schemes PUBLIC grid parameters malloc norms)
