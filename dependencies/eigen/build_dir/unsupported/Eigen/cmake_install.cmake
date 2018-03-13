# Install script for directory: /media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/AdolcForward"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/AlignedVector3"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/ArpackSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/AutoDiff"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/BVH"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/EulerAngles"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/FFT"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/IterativeSolvers"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/KroneckerProduct"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/MatrixFunctions"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/MoreVectorization"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/MPRealSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/NonLinearOptimization"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/NumericalDiff"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/OpenGLSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/Polynomials"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/Skyline"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/SparseExtra"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/SpecialFunctions"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/Splines"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

