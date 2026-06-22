cmake_minimum_required(VERSION 3.10)
project(SLAE_S)

add_definitions(-DNO_FONT_AWESOME)

# C++17
set(CMAKE_CXX_STANDARD 17)

# ШЛЯХИ до заголовків
include_directories(
    libs/eigen
)

# Список файлів
set(SOURCES
    main.cpp

    Vector.cpp
    Matrix.cpp

    Preconditioner.cpp
    Solvers.cpp

    Logger.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})
