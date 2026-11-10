file(GLOB PCP_UNIT_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/unit/*.cpp)
add_executable(pcp_tests test_main.cpp ${PCP_UNIT_TESTS})
target_link_libraries(pcp_tests PRIVATE pcp doctest::doctest pcp_options)
target_compile_definitions(pcp_tests PRIVATE PCP_MODELS_DIR="${PROJECT_SOURCE_DIR}/models")

add_test(NAME unit COMMAND pcp_tests --test-suite-exclude=models,profile)
set_tests_properties(unit PROPERTIES LABELS unit)

# End-to-end tests on the provided models (slow; run in Release).
foreach(name ispar clevis tailsitter tailsitter_beta)
  add_test(NAME model_${name} COMMAND pcp_tests -tc=model_${name})
  set_tests_properties(model_${name} PROPERTIES LABELS models TIMEOUT 7200)
endforeach()
add_test(NAME models_misc COMMAND pcp_tests -ts=models -tce=model_*)
set_tests_properties(models_misc PROPERTIES LABELS models)
