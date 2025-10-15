# User note:
# This aims to work almost identical to libtool-based compilation pipeline.
# At project root, run:
# $ mkdir build && cd build/
# $ cmake -DMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ..

# Minimum required CMake version
cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project/Package name and version
project(fplll VERSION 5.4.5)
set(FPLLL_MAJOR_VERSION ${fplll_VERSION_MAJOR})
set(FPLLL_MINOR_VERSION ${fplll_VERSION_MINOR})
set(FPLLL_MICRO_VERSION ${fplll_VERSION_PATCH})
set(FPLLL_VERSION ${fplll_VERSION})

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION MATCHES "^14.")
    message(STATUS "You're using ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  else()
    message(WARNING "!! We recommend that you use clang-14 for this bootcamp. You're using ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}, a different version.")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  message(STATUS "You're using ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
else()
  message(WARNING "!! We recommend that you use clang-14 for this bootcamp. You're using ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}, which is not clang.")
endif()

# standard library use c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_SYSROOT $ENV{SDKROOT})

# Find system dependencies
# GMP
find_path(GMP_INCLUDE_DIR gmp.h PATHS /opt/homebrew/include /opt/homebrew/opt/gmp/include)
find_library(GMP_LIBRARY NAMES gmp libgmp PATHS /opt/homebrew/lib /opt/homebrew/opt/gmp/lib)
if (GMP_INCLUDE_DIR AND GMP_LIBRARY)
  set(GMP_FOUND TRUE)
  set(HAVE_LIBGMP 1)
  set(LIBGMP_LIB "-lgmp")
  message(STATUS "Found GMP: ${GMP_LIBRARY}")
else ()
  set(GMP_FOUND FALSE)
  message(FATAL_ERROR "GMP not found")
endif ()
mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)

# MPFR
find_path(MPFR_INCLUDE_DIR mpfr.h PATHS /opt/homebrew/include /opt/homebrew/opt/mpfr/include)
find_library(MPFR_LIBRARY NAMES mpfr libmpfr PATHS /opt/homebrew/lib /opt/homebrew/opt/mpfr/lib)
if (MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
  set(MPFR_FOUND TRUE)
  set(HAVE_LIBMPFR 1)
  set(LIBMPFR_LIB "-lmpfr")
  message(STATUS "Found MPFR: ${MPFR_LIBRARY}")
else ()
  set(MPFR_FOUND FALSE)
  message(FATAL_ERROR "MPFR not found")
endif ()
mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

# QD (Optional)
find_path(QD_INCLUDE_DIR qd/dd_real.h PATHS /opt/homebrew/include /opt/homebrew/opt/qd/include)
find_library(QD_LIBRARY NAMES qd libqd PATHS /opt/homebrew/lib /opt/homebrew/opt/qd/lib)
if (QD_INCLUDE_DIR AND QD_LIBRARY)
  set(QD_FOUND TRUE)
  set(FPLLL_WITH_QD 1)
  set(LIBQD_LIB "-lqd")
  message(STATUS "Found QD: ${QD_LIBRARY}")
else ()
  set(QD_FOUND FALSE)
  message(WARNING "QD not found, proceeding without QD")
endif ()
mark_as_advanced(QD_INCLUDE_DIR QD_LIBRARY)

find_package(Threads REQUIRED)
include(CheckCXXSourceCompiles)
if(Threads_FOUND)
  set(HAVE_PTHREAD 1)
  set(PTHREAD_LIB "-lpthread")
  set(PTHREAD_CFLAGS "-pthread")
  check_cxx_source_compiles("
    #include <pthread.h>
    int main() {
        pthread_mutexattr_t attr;
        return pthread_mutexattr_setprotocol(&attr, PTHREAD_PRIO_INHERIT);
    }
  " HAVE_PTHREAD_PRIO_INHERIT)
else()
  message(STATUS "pthread headers and libraries not found")
endif()

# Define project specific variables
# TODO: (alex) revert
# set(FPLLL_MAX_ENUM_DIM 256)
set(FPLLL_MAX_ENUM_DIM 0)
set(FPLLL_WITH_RECURSIVE_ENUM 1)
# TODO: (alex) revert
# set(FPLLL_MAX_PARALLEL_ENUM_DIM 120)
set(FPLLL_MAX_PARALLEL_ENUM_DIM 0)

# Use external enumeration
option(WITH_EXTENUM_DIR "Specify external enumeration library install directory" OFF)
if(WITH_EXTENUM_DIR)
    link_directories(${WITH_EXTENUM_DIR})
endif()

option(WITH_EXTENUM_LIB "Specify external enumeration library name" OFF)
if(WITH_EXTENUM_LIB)
    set(EXTENUM_LIB ${WITH_EXTENUM_LIB})
endif()

option(WITH_EXTENUM_FUNC "Specify external enumeration function name" OFF)
if(WITH_EXTENUM_FUNC)
    set(EXTENUM_FUNC ${WITH_EXTENUM_FUNC})
endif()

if(EXTENUM_LIB)
    if(NOT EXTENUM_FUNC)
        set(EXTENUM_FUNC ${EXTENUM_LIB})
    endif()
    try_compile(EXTENUM_LINK_SUCCESS ${CMAKE_BINARY_DIR} SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/extenum_check.cpp
        LINK_LIBRARIES ${EXTENUM_LIB}
    )
    if(NOT EXTENUM_LINK_SUCCESS)
        message(FATAL_ERROR "Failed to find external enumeration")
    endif()
    set(EXTENUM_LIBS "-l${EXTENUM_LIB}")
    add_definitions(-DFPLLL_EXTENUM_FUNC=${EXTENUM_FUNC})
endif()
set(EXTENUM_LIBS "${EXTENUM_LIBS}" CACHE STRING "External enumeration libraries")

# Include directories
include_directories(${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR} ${QD_INCLUDE_DIR})
set(LIBS "${LIBGMP_LIB} ${LIBMPFR_LIB} ${LIBQD_LIB} ${PTHREAD_LIB}")
message(STATUS "LIBS: ${LIBS}")

# setup fplll/fplll_config.h
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/fplll/fplll_config.cmake.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/fplll/fplll_config.h
  @ONLY)

# setup ./config.h
include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckTypeSize)
set(CMAKE_REQUIRED_FLAGS "-isysroot ${CMAKE_SYSROOT}")
check_include_file_cxx(cstdio HAVE_CSTDIO)
check_include_file_cxx(iostream HAVE_IOSTREAM)
check_include_file_cxx(string HAVE_STRING)
check_include_file_cxx(vector HAVE_VECTOR)
check_include_file_cxx(cstdint HAVE_CXX11)
check_include_file_cxx(limits HAVE_LIMITS)
check_include_file_cxx(limits.h HAVE_LIMITS_H)
check_include_file_cxx(dlfcn.h HAVE_DLFCN_H)
check_include_file_cxx(inttypes.h HAVE_INTTYPES_H)
check_include_file_cxx(stdbool.h HAVE_STDBOOL_H)
check_include_file_cxx(stdint.h HAVE_STDINT_H)
check_include_file_cxx(stdio.h HAVE_STDIO_H)
check_include_file_cxx(stdlib.h HAVE_STDLIB_H)
check_include_file_cxx(strings.h HAVE_STRINGS_H)
check_include_file_cxx(string.h HAVE_STRING_H)
check_include_file_cxx(sys/stat.h HAVE_SYS_STAT_H)
check_include_file_cxx(sys/types.h HAVE_SYS_TYPES_H)
check_include_file_cxx(unistd.h HAVE_UNISTD_H)
check_function_exists(floor HAVE_FLOOR)
check_function_exists(pow HAVE_POW)
check_function_exists(rint HAVE_RINT)
check_function_exists(sqrt HAVE_SQRT)
check_function_exists(strtol HAVE_STRTOL)
check_type_size("_Bool" HAVE__BOOL LANGUAGE C)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h
  @ONLY)

# Set up fplll package, corresponding to fplll/Makefile.am and ./fplll.pc.in
# NOTE: ignored distribution related stuff, EXTRA_DIST not translated over
include(GNUInstallDirs)

set(FPLLL_PKG_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/fplll.pc")
# setup fplll.pc for package
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/fplll.cmake.pc.in
  ${FPLLL_PKG_CONFIG}
  @ONLY)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")
## enabled via -DDEBUG
if(DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# $(srcdir) in Makefile.am
set(FPLLL_SRCDIR "${CMAKE_CURRENT_SOURCE_DIR}/fplll")
# $(include_fpllldir) in Makefile.am
set(FPLLL_PKG_INCLUDEDIR "${CMAKE_INSTALL_FULL_INCLUDEDIR}/fplll")
# $(pkgconfigdir) in Makefile.am
set(FPLLL_PKG_CONFIGDIR "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig")
# $(pkgdatadir) in Makefile.am
set(FPLLL_PKG_DATADIR "${CMAKE_INSTALL_FULL_DATADIR}/fplll")
# $(strategydir) in Makefile.am
set(FPLLL_PKG_STRATEGYDIR "${FPLLL_PKG_DATADIR}/strategies")

# === libfplll, fplll headers ===
# FPLLL_LT_CURRENT in autoconf
set(LIBFPLLL_MAJOR_VERSION 8)
# FPLLL_LT_REVISION in autoconf
set(LIBFPLLL_MINOR_VERSION 1)
# FPLLL_LT_AGE in autoconf
set(LIBFPLLL_MICRO_VERSION 0)
# NOTE: this is library version, different from overall package/project version
set(LIBFPLLL_VERSION "${LIBFPLLL_MAJOR_VERSION}.${LIBFPLLL_MINOR_VERSION}.${LIBFPLLL_MICRO_VERSION}")

# Specify the additional files to be cleaned in the current directory
set(COVERAGE_ARTIFACTS "*.gcov" ".libs/*.gcda" ".libs/*.gcno" "*.gcno" "*.gcda")
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${COVERAGE_ARTIFACTS}")

set(LIBFPLLL_SOURCES
  fplll.cpp fplll.h
  util.cpp util.h
  enum/topenum.cpp enum/topenum.h
  enum/enumerate.cpp enum/enumerate.h
  enum/enumerate_base.cpp enum/enumerate_base.h
  enum/enumerate_ext.cpp enum/enumerate_ext.h enum/enumerate_ext_api.h
  enum/evaluator.cpp enum/evaluator.h
  lll.cpp lll.h
  wrapper.cpp wrapper.h
  bkz.cpp bkz.h
  bkz_param.cpp bkz_param.h
  gso_interface.cpp gso_interface.h gso_gram.cpp gso_gram.h gso.cpp gso.h
  pruner/pruner.cpp
  pruner/pruner.h
  pruner/pruner_simplex.h
  pruner/pruner_cost.cpp
  pruner/pruner_optimize.cpp
  pruner/pruner_optimize_tc.cpp
  pruner/pruner_optimize_tp.cpp
  pruner/pruner_prob.cpp
  pruner/pruner_util.cpp
  householder.cpp householder.h hlll.cpp hlll.h
  io/json.hpp
  threadpool.h threadpool.cpp io/thread_pool.hpp
)
if (FPLLL_MAX_PARALLEL_ENUM_DIM GREATER 0)
  list(APPEND LIBFPLLL_SOURCES
    enum-parallel/enumeration.h enum-parallel/enumlib.h enum-parallel/fplll_types.h
    enum-parallel/enumlib.cpp
    enum-parallel/enumlib_dim.160.cpp enum-parallel/enumlib_dim.150.cpp enum-parallel/enumlib_dim.140.cpp enum-parallel/enumlib_dim.130.cpp
    enum-parallel/enumlib_dim.120.cpp enum-parallel/enumlib_dim.110.cpp enum-parallel/enumlib_dim.100.cpp enum-parallel/enumlib_dim.90.cpp
    enum-parallel/enumlib_dim.80.cpp enum-parallel/enumlib_dim.70.cpp enum-parallel/enumlib_dim.60.cpp enum-parallel/enumlib_dim.50.cpp
    enum-parallel/enumlib_dim.40.cpp enum-parallel/enumlib_dim.30.cpp enum-parallel/enumlib_dim.20.cpp
  )
endif()
list(TRANSFORM LIBFPLLL_SOURCES PREPEND "fplll/")
# NOTE: EXTRA_libfplll_la_SOURCES+ are already included in the above sources
# e.g. svpcvp.cpp is included in fplll.cpp
# e.g. enumlib_dim.cpp is included in any enumlib_dim.x.cpp


# First define the (versioned) dynamic library
add_library(libfplll_dynlib SHARED ${LIBFPLLL_SOURCES})
add_library(libfplll_static STATIC ${LIBFPLLL_SOURCES})
if(DEBUG)
  add_library(libfplll_debug SHARED ${LIBFPLLL_SOURCES})
endif()

# libfplll_la_LIBADD in fplll/Makefile.am
set(LIBFPLLL_LIBS "${LIBS} ${EXTENUM_LIB}")
string(STRIP ${LIBFPLLL_LIBS} LIBFPLLL_LIBS) # remove leading and trailing whitespaces
target_link_libraries(libfplll_dynlib PRIVATE ${LIBFPLLL_LIBS})
target_link_libraries(libfplll_static PRIVATE ${LIBFPLLL_LIBS})
if(DEBUG)
  target_link_libraries(libfplll_debug PRIVATE ${LIBFPLLL_LIBS})
endif()

# libfplll_la_LDFLAGS in fplll/Makefile.am
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_link_options(libfplll_dynlib PRIVATE
    -Wl,-no-undefined
    -Wl,-version-info ${LIBFPLLL_VERSION}
    ${PTHREAD_CFLAGS}
  )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_link_options(libfplll_dynlib PRIVATE
    # MacOS's ld already err on undefined by default
    -Wl,-current_version,${LIBFPLLL_VERSION}
    ${PTHREAD_CFLAGS}
  )
endif()

# libfplll_la_CXXFLAGS in fplll/Makefile.am
separate_arguments(CXXFLAGS_LIST UNIX_COMMAND "${CMAKE_CXX_FLAGS} ${PTHREAD_CFLAGS}")
target_compile_options(libfplll_dynlib PRIVATE ${CXXFLAGS_LIST})
target_compile_options(libfplll_static PRIVATE ${CXXFLAGS_LIST})
if(DEBUG)
  separate_arguments(CXXFLAGS_DEBUG_LIST UNIX_COMMAND "-DDEBUG ${CMAKE_CXX_FLAGS} ${PTHREAD_CFLAGS}")
  target_compile_options(libfplll_debug PRIVATE ${CXXFLAGS_DEBUG_LIST})
endif()

set_target_properties(libfplll_dynlib PROPERTIES
  PREFIX "" # already have lib prefix in the name
  OUTPUT_NAME "libfplll"
  VERSION ${LIBFPLLL_VERSION}
)
set_target_properties(libfplll_static PROPERTIES
  PREFIX ""
  OUTPUT_NAME "libfplll"
  VERSION ${LIBFPLLL_VERSION}
)
if(DEBUG)
  set_target_properties(libfplll_debug PROPERTIES
    PREFIX ""
    OUTPUT_NAME "libfpllld"
    VERSION ${LIBFPLLL_VERSION}
  )
endif()


# header files to include/
install(
  DIRECTORY ${FPLLL_SRCDIR}
  TYPE INCLUDE
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "nr/matrix.cpp"
  PATTERN "*.hpp"
  PATTERN "*.inl"
  PATTERN "main.h" EXCLUDE
  PATTERN ".deps" EXCLUDE
  PATTERN ".libs" EXCLUDE
)
# dylib to lib/
install(TARGETS libfplll_dynlib libfplll_static DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR})

# fplll bin
add_executable(fplll fplll/main.cpp fplll/main.h config.h)
add_dependencies(fplll libfplll_static)
set(FPLLL_LIBS "${LIBS} ${EXTENUM_LIB} ${CMAKE_CURRENT_BINARY_DIR}/libfplll.a")
string(STRIP ${FPLLL_LIBS} FPLLL_LIBS) # remove leading and trailing whitespaces
target_link_libraries(fplll PRIVATE ${FPLLL_LIBS})
install(TARGETS fplll DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

if(DEBUG)
  add_executable(fplll_debug fplll/main.cpp fplll/main.h config.h)
  add_dependencies(fplll_debug libfplll_static)
  target_compile_options(fplll_debug PRIVATE ${CXXFLAGS_DEBUG_LIST})
  target_link_libraries(fplll_debug PRIVATE ${FPLLL_LIBS})
  install(TARGETS fplll_debug DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
endif()

# latticegen bin
add_executable(latticegen fplll/latticegen.cpp config.h)
add_dependencies(latticegen libfplll_static)
target_link_libraries(latticegen PRIVATE ${FPLLL_LIBS})
install(TARGETS latticegen DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

if(DEBUG)
  add_executable(latticegen_debug fplll/latticegen.cpp config.h)
  add_dependencies(latticegen_debug libfplll_static)
  target_compile_options(latticegen_debug PRIVATE ${CXXFLAGS_DEBUG_LIST})
  target_link_libraries(latticegen_debug PRIVATE ${FPLLL_LIBS})
  install(TARGETS latticegen_debug DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
endif()


# llldiff bin
add_executable(llldiff fplll/llldiff.cpp config.h)
add_dependencies(llldiff libfplll_static)
target_link_libraries(llldiff PRIVATE ${FPLLL_LIBS})
install(TARGETS llldiff DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

# Define preprocessor definitions
add_compile_definitions(
  FPLLL_DEFAULT_STRATEGY_PATH="${FPLLL_PKG_STRATEGYDIR}"
  FPLLL_DEFAULT_STRATEGY="${FPLLL_PKG_STRATEGYDIR}/default.json"
)
install(FILES ${FPLLL_PKG_CONFIG} DESTINATION ${FPLLL_PKG_CONFIGDIR})
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/strategies/default.json" DESTINATION ${FPLLL_PKG_STRATEGYDIR})

# tests
# ${TOPSRCDIR} in tests/Makefile.am
set(TOPSRCDIR ${CMAKE_CURRENT_SOURCE_DIR})
find_path(INSTALLED_FPLLL_H fplll.h PATHS "${FPLLL_PKG_INCLUDEDIR}")
if (INSTALLED_FPLLL_H)
  add_subdirectory(tests)
else()
  message(WARNING "run `make install` to install header files, then run `make` again")
endif()

# extra Makefile target
set(CLANGFORMAT "clang-format")
add_custom_target(check-style
    COMMAND ${CMAKE_COMMAND} -E echo "Running style checks..."
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/.check-m4.sh
    COMMAND ${CLANGFORMAT} -i --style=file ${CMAKE_CURRENT_SOURCE_DIR}/fplll/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/fplll/*.h ${CMAKE_CURRENT_SOURCE_DIR}/fplll/*/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/fplll/*/*.h ${CMAKE_CURRENT_SOURCE_DIR}/fplll/*/*.inl ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Checking code style with clang-format"
)
