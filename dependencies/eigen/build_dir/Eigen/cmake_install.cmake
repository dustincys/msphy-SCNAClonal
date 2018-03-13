# Install script for directory: /media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Cholesky"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/CholmodSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Core"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Dense"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Eigen"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Eigenvalues"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Geometry"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Householder"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/IterativeLinearSolvers"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Jacobi"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/LU"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/MetisSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/OrderingMethods"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/PardisoSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/PaStiXSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/QR"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/QtAlignedMalloc"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/Sparse"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SparseCholesky"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SparseCore"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SparseLU"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SparseQR"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SPQRSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/StdDeque"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/StdList"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/StdVector"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SuperLUSupport"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/SVD"
    "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/UmfPackSupport"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/media/d/github/phy-SCNAClonal/dependencies/eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

