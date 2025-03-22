
file(GLOB TEST_SOURCES "*.cc")

foreach(TEST_SOURCE ${TEST_SOURCES})
   get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
   add_executable(${TEST_NAME} ${TEST_SOURCE})
   target_link_libraries(${TEST_NAME}
      Honeycomb2
   )
   set_target_properties(${TEST_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/exe-tests"
   )

   # Add test
   add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
