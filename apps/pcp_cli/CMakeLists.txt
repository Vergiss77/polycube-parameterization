add_executable(pcp_cli main.cpp)
target_link_libraries(pcp_cli PRIVATE pcp CLI11::CLI11 pcp_options)
