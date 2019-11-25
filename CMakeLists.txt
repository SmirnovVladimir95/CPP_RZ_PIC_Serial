cmake_minimum_required(VERSION 3.15)
project(CPP_RZ_PIC)

set(CMAKE_CXX_STANDARD 17)

add_executable(CPP_RZ_PIC main.cpp Tools/Matrix.cpp Tools/Matrix.h Field/PoissonSolver.cpp Field/PoissonSolver.h Field/Field.cpp Field/Field.h Particles/Particles.cpp Particles/Particles.h Grid/Grid.cpp Grid/Grid.h Particles/Pusher.cpp Particles/Pusher.h Particles/Interpolation.cpp Particles/Interpolation.h)