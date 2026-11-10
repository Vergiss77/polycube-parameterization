# All third-party dependencies are fetched at configure time with pinned versions.
include(FetchContent)
set(FETCHCONTENT_QUIET ON)
# Some dependencies declare very old cmake_minimum_required versions, which CMake 4 rejects.
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()

# libigl (header-only mode) + its pinned Eigen + TetGen (copyleft module).
set(LIBIGL_USE_STATIC_LIBRARY OFF CACHE BOOL "" FORCE)
set(LIBIGL_COPYLEFT_TETGEN ON CACHE BOOL "" FORCE)
foreach(opt LIBIGL_BUILD_TESTS LIBIGL_BUILD_TUTORIALS LIBIGL_INSTALL LIBIGL_EMBREE LIBIGL_GLFW
            LIBIGL_IMGUI LIBIGL_OPENGL LIBIGL_PNG LIBIGL_PREDICATES LIBIGL_SPECTRA LIBIGL_XML
            LIBIGL_COPYLEFT_CGAL LIBIGL_COPYLEFT_COMISO LIBIGL_COPYLEFT_CORE
            LIBIGL_RESTRICTED_MATLAB LIBIGL_RESTRICTED_MOSEK LIBIGL_RESTRICTED_TRIANGLE)
  set(${opt} OFF CACHE BOOL "" FORCE)
endforeach()
FetchContent_Declare(libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG v2.6.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(libigl)

FetchContent_Declare(nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
  URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
set(JSON_BuildTests OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(nlohmann_json)

FetchContent_Declare(cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.2
  GIT_SHALLOW TRUE)
set(CLI11_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CLI11_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CLI11_BUILD_DOCS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(cli11)

if(PCP_BUILD_TESTS)
  FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE)
  set(DOCTEST_WITH_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(doctest)
endif()

if(PCP_BUILD_VIEWER)
  FetchContent_Declare(polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG v2.3.0
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(polyscope)
endif()
