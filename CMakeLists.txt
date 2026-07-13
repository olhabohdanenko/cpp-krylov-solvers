cmake_minimum_required(VERSION 3.16)
project(SLAE_S LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

execute_process(COMMAND spack location -i trilinos OUTPUT_VARIABLE TRILINOS_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND spack location -i petsc OUTPUT_VARIABLE PETSC_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND spack location -i eigen OUTPUT_VARIABLE EIGEN_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

find_package(Trilinos REQUIRED HINTS ${TRILINOS_DIR})
find_package(Eigen3 REQUIRED HINTS ${EIGEN_DIR})
find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)

set(PRECOND_SOURCES
    src/Preconditioners/CustomPreconditioner.cpp
    src/Preconditioners/PetscPreconditioner.cpp
    src/Preconditioners/TrilinosPreconditioner.cpp
)

set(SOLVERS_SOURCES
    src/Solvers/CustomSolver.cpp
    src/Solvers/PetscSolver.cpp
    src/Solvers/TrilinosSolver.cpp
)

add_executable(SLAE_S 
    main.cpp
    src/Vector.cpp
    src/Matrix.cpp
    src/Logger.cpp
    src/MemoryTracker.cpp
    ${PRECOND_SOURCES}
    ${SOLVERS_SOURCES}
)

target_include_directories(SLAE_S PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_PATH}
)

target_link_libraries(SLAE_S PRIVATE
    PkgConfig::PETSC
    Trilinos::all_libs
    Kokkos::kokkos
    KokkosKernels::kokkoskernels
    Teuchos::all_libs
    MPI::MPI_CXX
)

target_link_directories(SLAE_S PRIVATE 
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
)

message(STATUS "=== Trilinos Debug Info ===")
message(STATUS "Trilinos_LIBRARIES: ${Trilinos_LIBRARIES}")
message(STATUS "Trilinos_TPL_LIBRARIES: ${Trilinos_TPL_LIBRARIES}")
message(STATUS "Trilinos_INCLUDE_DIRS: ${Trilinos_INCLUDE_DIRS}")
message(STATUS "Tpetra_LIBRARIES: ${Tpetra_LIBRARIES}")
message(STATUS "Kokkos_DIR: ${Kokkos_DIR}")
