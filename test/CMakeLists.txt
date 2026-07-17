# Unit testing
set(
  tests
  "tblite"
  "gfnff"
  "gfn0"
  "gfn0occ"
  "CN"
  "optimization"
)
set(
  test-srcs
  "testmol.f90"
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.F90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "lib-${PROJECT_NAME}-static"
  "test-drive::test-drive"
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()


# Focused COM-CV test target. It remains separate from CREST 3.0.1's legacy
# aggregate test driver, whose unrelated get_testmol_ link failure prevents a
# full aggregate link on this toolchain. Existing suites remain unchanged.
add_executable(
  "${PROJECT_NAME}-comcv-tester"
  "main_comcv.f90"
  "test_comcv.F90"
)
target_link_libraries(
  "${PROJECT_NAME}-comcv-tester"
  PRIVATE
  "lib-${PROJECT_NAME}-static"
  "test-drive::test-drive"
)
add_test(
  "${PROJECT_NAME}/comcv"
  "${PROJECT_NAME}-comcv-tester"
  "comcv"
)
