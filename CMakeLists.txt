cmake_minimum_required(VERSION 3.28)
project(N_Body)

set(CMAKE_CXX_STANDARD 23)

add_executable(N_Body main.cpp
        src/particle.h
        src/simulation.h)
