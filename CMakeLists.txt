cmake_minimum_required(VERSION 3.13)

project(openFEM
  DESCRIPTION "A free and open source quasi-static FEM extractor for circuit parameters (R/L/C)"
  LANGUAGES C)

# MSVC は C99 の一部 (VLA) を実装しないため C11 を指定する。
# ソースは C99 相当の範囲で書き、VLA は使わない (OpenFDTD と同じ規則)。
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ── Windows (MSVC) 対応 ──────────────────────────────────────────────
# libm は存在しない (MATH_LIB を空にする)。ソースは UTF-8 (日本語コメント)。
if(MSVC)
    add_compile_options(/utf-8)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS NOMINMAX)
    set(MATH_LIB "")
else()
    set(MATH_LIB m)
endif()

option(WITH_OPENMP "Enable OpenMP support" ON)
if(WITH_OPENMP)
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        message(STATUS "OpenMP enabled")
    else()
        message(STATUS "OpenMP not found (building single thread)")
    endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
foreach(cfg IN ITEMS DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${cfg} ${CMAKE_SOURCE_DIR}/bin)
endforeach()

include_directories(${CMAKE_SOURCE_DIR}/include)

set(SOL_SOURCES
    src/sol_Main.c
    sol/input_data.c
    sol/setup.c
    sol/ingeometry.c
    sol/crs.c
    sol/assemble.c
    sol/solver_cg.c
    sol/solver_cocg.c
    sol/solve.c
    sol/matutil.c
    sol/outputRLC.c
    sol/utils.c
)

add_executable(ofe ${SOL_SOURCES})
target_link_libraries(ofe PRIVATE ${MATH_LIB})
if(WITH_OPENMP AND OpenMP_C_FOUND)
    target_link_libraries(ofe PRIVATE OpenMP::OpenMP_C)
endif()

add_executable(ofe_post post/post_Main.c)
target_link_libraries(ofe_post PRIVATE ${MATH_LIB})

# 深い再帰は無いが、OpenFDTD と同じく MSVC のスタックを広げておく
if(MSVC)
    set_target_properties(ofe      PROPERTIES LINK_FLAGS "/STACK:16777216")
    set_target_properties(ofe_post PROPERTIES LINK_FLAGS "/STACK:16777216")
endif()

install(TARGETS ofe ofe_post RUNTIME DESTINATION bin)
