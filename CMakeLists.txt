cmake_minimum_required(VERSION 3.16)

project(DynamicMatrices LANGUAGES CXX)

set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)
set (PATCH_LEVEL 0)

# require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)




#the code subdir
add_subdirectory(src)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/DynamicMatricesConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DynamicMatrices
)

set(DynamicMatrices_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${PATCH_LEVEL})

#set_property(TARGET DynamicMatrices PROPERTY VERSION ${DynamicMatrices_VERSION})
#set_property(TARGET DynamicMatrices PROPERTY SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR})
set_property(TARGET DynamicMatrices PROPERTY INTERFACE_DynamicMatrices_MAJOR_VERSION ${VERSION_MAJOR})
set_property(TARGET DynamicMatrices APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING DynamicMatrices_MAJOR_VERSION
)

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/DynamicMatricesConfigVersion.cmake"
  VERSION "${DynamicMatrices_VERSION}"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/DynamicMatricesConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/DynamicMatricesConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DynamicMatrices
)

set(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/DynamicMatrices CACHE PATH "Installation directory for CMake files")

install(EXPORT DynamicMatricesTargets
        FILE DynamicMatricesTargets.cmake
        NAMESPACE DynamicMatrices::
        DESTINATION ${INSTALL_CMAKE_DIR}
)

#the tests subdir
add_subdirectory(tests)