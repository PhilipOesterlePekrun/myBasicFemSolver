cmake_minimum_required(VERSION 3.25)
project(myBasicFemSolver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Source Collection ---
file(GLOB_RECURSE SRC_FILES "src/*.cpp" "src/*.hpp")
file(GLOB_RECURSE TEST_FILES "test/*.cpp" "test/*.hpp")
list(FILTER SRC_FILES EXCLUDE REGEX ".*/_EXCLUDE_.*")
list(FILTER TEST_FILES EXCLUDE REGEX ".*/_EXCLUDE_.*")



# --- Self-include ---
include_directories(
  "src"
  "src/GlobalDefinitions"
)

# --- SFML ---
get_filename_component(SFML_ROOT "$ENV{HOME}/misc/tpl/SFML_base/install" ABSOLUTE)
set(CMAKE_PREFIX_PATH "${SFML_ROOT}" CACHE PATH "SFML install prefix" FORCE)
include_directories("${SFML_ROOT}/include")
link_directories("${SFML_ROOT}/lib")
find_package(SFML REQUIRED COMPONENTS Graphics Window System)

# --- Add Executable ---
add_executable(myBasicFemSolver ${SRC_FILES})

# --- Create executable for testing/demo ---
add_executable(UserTest test/user_test/UserTest.cpp)

# --- myUtils ---
set(MYUTILS_DIR "$ENV{HOME}/misc/myUtils_Base/install")
option(myUtils_ENABLE_TIMERS "Enable timers" ON)
add_compile_definitions(
  myUtils_ENABLE_TIMERS=$<IF:$<BOOL:${myUtils_ENABLE_TIMERS}>,1,0>
)
option(myUtils_DbPr_ON "Enable timers" ON)
add_compile_definitions(
  myUtils_DbPr_ON=$<IF:$<BOOL:${myUtils_DbPr_ON}>,1,0>
)

target_include_directories(myBasicFemSolver PRIVATE ${MYUTILS_DIR}/include)
target_include_directories(UserTest PRIVATE ${MYUTILS_DIR}/include "src")





target_link_libraries(myBasicFemSolver
    PRIVATE
        ${MYUTILS_DIR}/lib/libmyUtils.a
        sfml-graphics sfml-window sfml-system
)
target_link_libraries(UserTest
    PRIVATE
        ${MYUTILS_DIR}/lib/libmyUtils.a
        sfml-graphics sfml-window sfml-system
)
