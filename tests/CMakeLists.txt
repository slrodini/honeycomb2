
# file(GLOB TEST_SOURCES "*.cc")

configure_file(
   ${CMAKE_CURRENT_SOURCE_DIR}/test_config.hpp.in
   ${CMAKE_CURRENT_BINARY_DIR}/test_config.hpp
   @ONLY
)

set(TEST_SOURCES  test_alpha_s.cc 
                  test_S_to_T.cc
                  test_evolution_operator.cc
                  test_foreign_interface.cc
                  # test_complete_evolution.cc
)

foreach(TEST_SOURCE ${TEST_SOURCES})
   get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
   add_executable(${TEST_NAME} ${TEST_SOURCE})
   target_include_directories(${TEST_NAME} PRIVATE
      ${PROJECT_BINARY_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/inc
   )

   target_link_libraries(${TEST_NAME}
      Honeycomb2
   )
   set_target_properties(${TEST_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/exe-tests"
   )

   # Add test
   # add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
   message(STATUS "Adding test for ${TEST_NAME}")
   add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})

endforeach()
