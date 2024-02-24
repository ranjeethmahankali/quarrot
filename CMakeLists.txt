cmake_minimum_required(VERSION 3.20.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# vcpkg
if(NOT DEFINED ENV{VCPKG_PATH})
  message(FATAL_ERROR "Please clone vcpkg and set the VCPKG_PATH environment variable")
endif()

set(VCPKG_FEATURE_FLAGS versions)
set(VCPKG_TARGET_TRIPLET x64-linux)
set(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_PATH}/scripts/buildsystems/vcpkg.cmake)

project(quarrot)

find_package(OpenMesh REQUIRED)
find_package(glm CONFIG REQUIRED)
add_compile_definitions(GLM_ENABLE_EXPERIMENTAL)

file(GLOB SRC_FILES "src/*.cpp")
add_library(quarrot STATIC ${SRC_FILES})
target_link_libraries(quarrot PUBLIC OpenMeshCore OpenMeshTools glm::glm)
target_include_directories(quarrot PUBLIC "include/")

# Testing
find_package(Catch2 CONFIG REQUIRED)
file(GLOB TEST_FILES "tests/*.cpp")
add_executable(test_quarrot ${TEST_FILES})
target_include_directories(test_quarrot PRIVATE "include/")
target_link_libraries(test_quarrot PRIVATE quarrot Catch2::Catch2 Catch2::Catch2WithMain)
include(CTest)
include(Catch)
catch_discover_tests(test_quarrot)
