cmake_minimum_required(VERSION 3.16)

project(clusteringPPTests LANGUAGES CXX)

# require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# enable FetchContent to get Catch2
include(FetchContent)

# declare Catch2
FetchContent_Declare(Catch2                              # name of the content
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git  # the repository
  GIT_TAG        v3.0.1                               # the tag
  )

# make available
FetchContent_MakeAvailable(Catch2)

#this should make the tests standalone
if(NOT TARGET DynamicMatrices::DynamicMatrices)
  MESSAGE("Not found DynamicMatrices::DynamicMatrices")
  #find_package(DynamicMatrices REQUIRED)
endif()

# testing binary
add_executable(DM_test MatricesTests.cpp)
target_link_libraries(DM_test PRIVATE DynamicMatrices::DynamicMatrices Catch2::Catch2WithMain )
target_include_directories(DM_test
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )



# enable testing functionality
enable_testing()

include(CTest)
message(${catch2_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(DM_test)