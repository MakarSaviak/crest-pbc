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
  "helpers.f90"
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

# This validation-only target exercises the process-MTD capsule format,
# including explicit rejection of a V4 header by the V5 reader.  It is a CTest
# regression and is not part of normal production MTD execution.
if(WITH_GFNFF)
  add_executable(
    "${PROJECT_NAME}-mtd-capsule-regression"
    "${PROJECT_SOURCE_DIR}/prototype_tests/capsule_roundtrip_driver.F90"
  )
  target_link_libraries(
    "${PROJECT_NAME}-mtd-capsule-regression"
    PRIVATE
    "lib-${PROJECT_NAME}-static"
  )
  add_test(
    "${PROJECT_NAME}/mtd-capsule-v5"
    "${PROJECT_NAME}-mtd-capsule-regression"
    "${CMAKE_CURRENT_BINARY_DIR}/mtd-capsule-v5.bin"
  )
endif()
