cmake_minimum_required(VERSION 3.16)
project(atom)

set(CMAKE_RELEASE_CXX_FLAGS -O3)
set(CMAKE_DEBUG_CXX_FLAGS -O0)

if( UNIX AND NOT APPLE)
        set(CMAKE_CXX_STANDARD 23)
else()
        set(CMAKE_CXX_STANDARD 20)
endif()


include_directories(
        src/
        src/stdeig
        src/montecarlo
        src/montecarlo/mersenne
        src/montecarlo/mrg32k3a
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
if( UNIX AND NOT APPLE)
        target_link_libraries(twoCircles pthread)
endif()

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

add_executable(secondIntegral_a src/montecarlo/SecondIntegral_a.cpp ${MONTE_CARLO_LIB_SRC})
if( UNIX AND NOT APPLE)
        target_link_libraries(secondIntegral_a pthread)
endif()

add_executable(secondIntegral_a2 src/montecarlo/SecondIntegral_a2.cpp ${MONTE_CARLO_LIB_SRC})
if( UNIX AND NOT APPLE)
        target_link_libraries(secondIntegral_a2 pthread)
endif()

if( UNIX AND NOT APPLE)
        include_directories(
                src/montecarlo/cuda/cudalib
        )

        set(CUDA_SRC
                src/montecarlo/cuda/cudalib/CudaArray.h
                src/montecarlo/cuda/cudalib/CudaModule.h
                src/montecarlo/cuda/cudalib/CudaFunction.h
                src/montecarlo/cuda/cudalib/CudaDevice.cpp
                src/montecarlo/cuda/cudalib/CudaDevice.h
                src/montecarlo/cuda/cudalib/CudaContext.cpp
                src/montecarlo/cuda/cudalib/CudaContext.h src/montecarlo/cuda/FirstIntegralCuda.cpp)

        add_executable(basicCudaTest src/montecarlo/cuda/BasicCudaTest.cpp ${CUDA_SRC} src/montecarlo/mersenne/TestMersenneTwister.cpp src/montecarlo/mrg32k3a/TestMrg32k3a.cpp)
        target_link_libraries(basicCudaTest pthread cuda)

        add_executable(firstIntegralCuda src/montecarlo/cuda/FirstIntegralCuda.cpp ${CUDA_SRC})
        target_link_libraries(firstIntegralCuda pthread cuda)
endif()

set(MERSENNE_SRC
        src/montecarlo/mersenne/MersenneTwister.h
        src/montecarlo/mersenne/DenseF2Matrix.h
        src/montecarlo/mersenne/SparseMatrix.h
        src/montecarlo/mersenne/MatrixElement.h
        )
add_executable(testMersenneTwister src/montecarlo/mersenne/TestMersenneTwister.cpp ${MERSENNE_SRC})

add_executable(testJumpAhead src/montecarlo/mersenne/TestJumpAhead.cpp ${MERSENNE_SRC})

add_executable(testCharacteristicPoly src/montecarlo/mersenne/TestCharacteristicPoly.cpp ${MERSENNE_SRC})

set(MRG32K3A_SRC
        src/montecarlo/mrg32k3a/Matrix3x3.h
        src/montecarlo/mrg32k3a/Mrg32k3a.h
        )

add_executable(testMrg32k3a src/montecarlo/mrg32k3a/TestMrg32k3a.cpp ${MRG32K3A_SRC})
if( UNIX AND NOT APPLE)
        target_link_libraries(testMrg32k3a pthread)
endif()
