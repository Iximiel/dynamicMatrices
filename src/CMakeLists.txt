cmake_minimum_required(VERSION 3.16)

project(DynamicMatrices LANGUAGES CXX)

# require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF(${OPT_DEBUGWARNINGS})
    if(MSVC)
      add_compile_options(/W4 /WX)
    else()
      add_compile_options(
          -Wall
          -Wextra
          -Wshadow
          -Wnon-virtual-dtor
          -Wold-style-cast
          -Wcast-align
          -Wunused
          -Woverloaded-virtual
          -Wpedantic
          -Wconversion
          -Wsign-conversion
          -Wnull-dereference
          -Wdouble-promotion
          -Wformat=2
          -Wduplicated-cond
          -Wduplicated-branches
          -Wlogical-op
          -Wuseless-cast
          #-Wlifetime
          )
    endif()
endif()

add_library(DynamicMatrices INTERFACE)
add_library(DynamicMatrices::DynamicMatrices ALIAS DynamicMatrices)

target_include_directories(DynamicMatrices
  INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )

install(TARGETS DynamicMatrices
 EXPORT DynamicMatricesTargets
 RUNTIME DESTINATION ${INSTALL_BIN_DIR}
 LIBRARY DESTINATION ${INSTALL_LIB_DIR}
 ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
 PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
 INCLUDES DESTINATON ${INSTALL_INCLUDE_DIR}
 )


 