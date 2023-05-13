cmake_minimum_required (VERSION 2.8)
project (tinyraytracer)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
message(STATUS "SRC_DIR: ${SRC_DIR}")
message(STATUS "LIB_DIR: ${LIB_DIR}")

include(CheckCXXCompilerFlag)
function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-pedantic")
enable_cxx_compiler_flag_if_supported("-std=c++17")

enable_cxx_compiler_flag_if_supported("-O3")
enable_cxx_compiler_flag_if_supported("-fopenmp")


## Use OpenMP in macOS

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Clang detected")
    enable_cxx_compiler_flag_if_supported("-Xpreprocessor -fopenmp")
    enable_cxx_compiler_flag_if_supported("-lomp")
endif()
# find_package(OpenMP REQUIRED)
# target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
## or
# set(OPENMP_LIBRARIES "/opt/homebrew/opt/libomp/lib")
# set(OPENMP_INCLUDES "/opt/homebrew/opt/libomp/include")
# include_directories(${OPENMP_INCLUDES})
# link_directories(${OPENMP_LIBRARIES})


file(GLOB SOURCES
    "${SRC_DIR}/*.cpp"
)
# set(SOURCES
#     "${SRC_DIR}/main.cpp"
#     "${SRC_DIR}/model.cpp"
#     "${SRC_DIR}/utils.cpp"
#     )

# Executable definition and properties
add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE "${SRC_DIR}")
target_include_directories(${PROJECT_NAME} PRIVATE "${LIB_DIR}")

# LIB_DIR
set(LIB_DIR "${LIB_DIR}/lib")
target_include_directories(${PROJECT_NAME} PRIVATE "${LIB_DIR}")
