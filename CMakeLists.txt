cmake_minimum_required(VERSION 3.16)
project(atom)

set(CMAKE_RELEASE_CXX_FLAGS -O3)
set(CMAKE_DEBUG_CXX_FLAGS -O0)

set(CMAKE_CXX_STANDARD 20)

include_directories(
        src/
        src/stdeig
        src/montecarlo
        )

set(PHYSICS_SRC
        src/grid/RealScalarPhysics.h
        src/grid/DiracPhysics.h
        src/stdeig/EigenSolver.h
        src/grid/DiracPhysics.h
        src/stdeig/ChebyshevFilter.cpp
        src/stdeig/ChebyshevFilter.h
        src/Constants.h
        src/grid/RealScalarAntisymmetricPhysics.h
        src/Energy.h
        )

add_executable(hydrogen src/Hydrogen.cpp)
target_link_libraries(hydrogen blas lapacke)

add_executable(dirac src/Dirac.cpp)
target_link_libraries(dirac blas lapacke)

add_executable(singleparticle_fd
        src/manybody1D/singleparticle/SingleParticleFiniteDifference.cpp
        src/manybody1D/singleparticle/ForceCurve.h
        src/manybody1D/singleparticle/Hamiltonian.h
        src/manybody1D/singleparticle/Snapshot.h
        ${PHYSICS_SRC}
        )
target_link_libraries(singleparticle_fd blas lapacke)

add_executable(twoparticle_fd src/manybody1D/twoparticle/TwoParticleFiniteDifference.cpp ${PHYSICS_SRC})
target_link_libraries(twoparticle_fd blas lapacke)

add_executable(twoCircles src/montecarlo/TwoCircles.cpp)

set(MONTE_CARLO_LIB_SRC
        src/montecarlo/SphericalPoint.h
        src/montecarlo/Bounds.h
        src/montecarlo/Generator.h
        src/montecarlo/PartialWork.h
        src/montecarlo/CommandLine.h
        src/montecarlo/BigSum.h
        src/montecarlo/EuclideanPoint.h
        src/MathUtil.h)

add_executable(spheres src/montecarlo/Spheres.cpp)

add_executable(testBigSum src/montecarlo/TestBigSum.cpp)

add_executable(firstIntegral src/montecarlo/FirstIntegral.cpp ${MONTE_CARLO_LIB_SRC})
if( UNIX AND NOT APPLE)
        target_link_libraries(firstIntegral pthread)
endif()

add_executable(secondIntegral src/montecarlo/SecondIntegral.cpp ${MONTE_CARLO_LIB_SRC})
