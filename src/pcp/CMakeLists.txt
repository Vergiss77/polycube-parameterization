file(GLOB_RECURSE PCP_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE PCP_HEADERS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
add_library(pcp STATIC ${PCP_SOURCES} ${PCP_HEADERS})
target_include_directories(pcp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(pcp
  PUBLIC igl::core Eigen3::Eigen nlohmann_json::nlohmann_json
  PRIVATE igl_copyleft::tetgen pcp_options)
