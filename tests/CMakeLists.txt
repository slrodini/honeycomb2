find_package(GTest)

if(GTest_FOUND)
   file(GLOB TEST_SOURCES "*.cc")

   foreach(TEST_SOURCE ${TEST_SOURCES})
      get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

      add_executable(${TEST_NAME} ${TEST_SOURCE})

      target_link_libraries(${TEST_NAME}
         Honeycomb2
         GTest::GTest
         GTest::Main
      )

      # Add test
      add_test(${TEST_NAME} ${TEST_NAME})
   endforeach()
else()
   message(WARNING "GTest not found - tests will not be built")
endif()