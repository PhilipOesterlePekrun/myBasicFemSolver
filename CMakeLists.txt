cmake_minimum_required(VERSION 3.25)
project(myBasicFemSolver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Source Collection ---
file(GLOB_RECURSE ALL_SRC_FILES "src/*.cpp" "src/*.hpp")

set(SRC_FILES "")
foreach(file ${ALL_SRC_FILES})
    get_filename_component(fname ${file} NAME)
    if(NOT fname MATCHES "^EXCLUDE.*")
        list(APPEND SRC_FILES ${file})
    endif()
endforeach()



# --- Self-include ---
include_directories("src/GlobalDefinitions")

# --- SFML ---
get_filename_component(SFML_ROOT "$ENV{HOME}/misc/tpl/SFML_base/install_SFML" ABSOLUTE)
set(CMAKE_PREFIX_PATH "${SFML_ROOT}" CACHE PATH "SFML install prefix" FORCE)
include_directories("${SFML_ROOT}/include")
link_directories("${SFML_ROOT}/lib")
find_package(SFML REQUIRED COMPONENTS Graphics Window System)

# --- Add Executable ---
add_executable(myBasicFemSolver ${SRC_FILES})

# --- myUtils ---
set(MYUTILS_DIR "$ENV{HOME}/util/myUtils_Base/install_myUtils")
target_include_directories(myBasicFemSolver PRIVATE ${MYUTILS_DIR}/include)





target_link_libraries(myBasicFemSolver
    PRIVATE
        ${MYUTILS_DIR}/lib/libmyUtils.a
        sfml-graphics sfml-window sfml-system
)
