# Warning and platform flags shared by all project targets.
add_library(pcp_options INTERFACE)
if(MSVC)
  target_compile_options(pcp_options INTERFACE /W3 /bigobj /permissive- /utf-8 /EHsc)
  target_compile_definitions(pcp_options INTERFACE NOMINMAX _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(pcp_options INTERFACE -Wall -Wextra -Wno-unused-parameter)
endif()
