cmake_minimum_required (VERSION 2.6)

project (ToolFrameworkApplicaiton)

set(DATAMODEL_INC "")
set(DATAMODEL_LIB_PATH "")
set(DATAMODEL_LIBS "")

set(MYTOOLS_INC "")
set(MYTOOLS_LIB_PATH "")
set(MYTOOLS_LIBS "")

include_directories (${DATAMODEL_INC} ${MYTOOLS_INC})
link_directories(${DATAMODEL_LIB_PATH} ${MYTOOLS_LIB_PATH})

MACRO(HEADER_DIRECTORIES return_list)
    FILE(GLOB_RECURSE new_list ${PROJECT_SOURCE_DIR}/*.h)
    FILE(COPY ${new_list} DESTINATION ${PROJECT_BINARY_DIR}/include)
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

FILE(COPY ${PROJECT_SOURCE_DIR}/configfiles DESTINATION ${PROJECT_BINARY_DIR}/)

HEADER_DIRECTORIES(header_list)
include_directories(${header_list})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

file(GLOB_RECURSE STORE_SRC RELATIVE ${CMAKE_SOURCE_DIR} "Dependencies/ToolFrameworkCore/src/Store/*.cpp")
add_library(Store SHARED ${STORE_SRC})

file(GLOB_RECURSE LOGGING_SRC RELATIVE ${CMAKE_SOURCE_DIR} "Dependencies/ToolFrameworkCore/src/Logging/*.cpp")
add_library(Logging SHARED ${LOGGING_SRC})

file(GLOB_RECURSE DATAMODEL_SRC RELATIVE ${CMAKE_SOURCE_DIR} "DataModel/*.cpp")
add_library(DataModel SHARED ${DATAMODEL_SRC})

file(GLOB_RECURSE MYTOOLS_SRC RELATIVE ${CMAKE_SOURCE_DIR} "UserTools/*.cpp")
add_library(MyTools SHARED ${MYTOOLS_SRC})

file(GLOB_RECURSE TOOLCHAIN_SRC RELATIVE ${CMAKE_SOURCE_DIR} "Dependencies/ToolFrameworkCore/src/ToolChain/*.cpp")
add_library(ToolChain SHARED ${TOOLCHAIN_SRC})

add_executable (main ${PROJECT_SOURCE_DIR}/src/main.cpp)
target_link_libraries (main Store Logging ToolChain MyTools DataModel pthread ${DATAMODEL_LIBS} ${MYTOOLS_LIBS})
