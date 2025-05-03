enable_testing()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${COVERAGE_ARTIFACTS}")

# -I all ${prefix}/include/fplll headers
include_directories(${CMAKE_INSTALL_FULL_INCLUDEDIR})
add_compile_definitions(TESTDATADIR="${TOPSRCDIR}")

add_executable(test_nr test_nr.cpp)
target_link_libraries(test_nr PRIVATE ${FPLLL_LIBS})
add_test(NAME TestNR COMMAND ${CMAKE_BINARY_DIR}/test_nr)

add_executable(test_lll EXCLUDE_FROM_ALL test_lll.cpp)
target_link_libraries(test_lll PRIVATE ${FPLLL_LIBS})
add_test(NAME TestLLL COMMAND ${CMAKE_BINARY_DIR}/test_lll)

add_executable(test_enum EXCLUDE_FROM_ALL test_enum.cpp)
target_link_libraries(test_enum PRIVATE ${FPLLL_LIBS})
add_test(NAME TestENUM COMMAND ${CMAKE_BINARY_DIR}/test_enum)

add_executable(test_svp EXCLUDE_FROM_ALL test_svp.cpp)
target_link_libraries(test_svp PRIVATE ${FPLLL_LIBS})
add_test(NAME TestSVP COMMAND ${CMAKE_BINARY_DIR}/test_svp)

add_executable(test_cvp EXCLUDE_FROM_ALL test_cvp.cpp)
target_link_libraries(test_cvp PRIVATE ${FPLLL_LIBS})
add_test(NAME TestCVP COMMAND ${CMAKE_BINARY_DIR}/test_cvp)

add_executable(test_bkz EXCLUDE_FROM_ALL test_bkz.cpp)
target_link_libraries(test_bkz PRIVATE ${FPLLL_LIBS})
add_test(NAME TestBKZ COMMAND ${CMAKE_BINARY_DIR}/test_bkz)

add_executable(test_pruner EXCLUDE_FROM_ALL test_pruner.cpp)
target_link_libraries(test_pruner PRIVATE ${FPLLL_LIBS})
add_test(NAME TestPRUNER COMMAND ${CMAKE_BINARY_DIR}/test_pruner)

add_executable(test_gso EXCLUDE_FROM_ALL test_gso.cpp)
target_link_libraries(test_gso PRIVATE ${FPLLL_LIBS})
add_test(NAME TestGSO COMMAND ${CMAKE_BINARY_DIR}/test_gso)

add_executable(test_lll_gram EXCLUDE_FROM_ALL test_lll_gram.cpp)
target_link_libraries(test_lll_gram PRIVATE ${FPLLL_LIBS})
add_test(NAME TestLLL_GRAM COMMAND ${CMAKE_BINARY_DIR}/test_lll_gram)

add_executable(test_hlll EXCLUDE_FROM_ALL test_hlll.cpp)
target_link_libraries(test_hlll PRIVATE ${FPLLL_LIBS})
add_test(NAME TestHLLL COMMAND ${CMAKE_BINARY_DIR}/test_hlll)

add_executable(test_svp_gram EXCLUDE_FROM_ALL test_svp_gram.cpp)
target_link_libraries(test_svp_gram PRIVATE ${FPLLL_LIBS})
add_test(NAME TestSVP_GRAM COMMAND ${CMAKE_BINARY_DIR}/test_svp_gram)

add_executable(test_bkz_gram EXCLUDE_FROM_ALL test_bkz_gram.cpp)
target_link_libraries(test_bkz_gram PRIVATE ${FPLLL_LIBS})
add_test(NAME TestBKZ_GRAM COMMAND ${CMAKE_BINARY_DIR}/test_bkz_gram)

add_executable(test_counter EXCLUDE_FROM_ALL test_counter.cpp)
target_link_libraries(test_counter PRIVATE ${FPLLL_LIBS})
add_test(NAME TestCOUNTER COMMAND ${CMAKE_BINARY_DIR}/test_counter)

add_executable(test_babai EXCLUDE_FROM_ALL test_babai.cpp)
target_link_libraries(test_babai PRIVATE ${FPLLL_LIBS})
add_test(NAME TestBABAI COMMAND ${CMAKE_BINARY_DIR}/test_babai)

add_executable(test_ceil EXCLUDE_FROM_ALL test_ceil.cpp)
target_link_libraries(test_ceil PRIVATE ${FPLLL_LIBS})
add_test(NAME TestCEIL COMMAND ${CMAKE_BINARY_DIR}/test_ceil)


add_custom_target(build_tests
  ALL DEPENDS test_nr test_lll test_enum
  test_svp test_cvp test_bkz test_pruner test_gso
  test_lll_gram test_hlll test_svp_gram test_bkz_gram
  test_counter test_babai test_ceil
)

# Create a custom target to run the tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS build_tests
)

# 'make test' also run 'ctest'
add_custom_target(test DEPENDS run_tests)
