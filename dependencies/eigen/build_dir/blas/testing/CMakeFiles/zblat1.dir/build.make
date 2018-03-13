# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/d/github/phy-SCNAClonal/dependencies/eigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir

# Include any dependencies generated for this target.
include blas/testing/CMakeFiles/zblat1.dir/depend.make

# Include the progress variables for this target.
include blas/testing/CMakeFiles/zblat1.dir/progress.make

# Include the compile flags for this target's objects.
include blas/testing/CMakeFiles/zblat1.dir/flags.make

blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o: blas/testing/CMakeFiles/zblat1.dir/flags.make
blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o: ../blas/testing/zblat1.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /media/d/github/phy-SCNAClonal/dependencies/eigen/blas/testing/zblat1.f -o CMakeFiles/zblat1.dir/zblat1.f.o

blas/testing/CMakeFiles/zblat1.dir/zblat1.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/zblat1.dir/zblat1.f.i"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /media/d/github/phy-SCNAClonal/dependencies/eigen/blas/testing/zblat1.f > CMakeFiles/zblat1.dir/zblat1.f.i

blas/testing/CMakeFiles/zblat1.dir/zblat1.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/zblat1.dir/zblat1.f.s"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /media/d/github/phy-SCNAClonal/dependencies/eigen/blas/testing/zblat1.f -o CMakeFiles/zblat1.dir/zblat1.f.s

blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.requires:

.PHONY : blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.requires

blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.provides: blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.requires
	$(MAKE) -f blas/testing/CMakeFiles/zblat1.dir/build.make blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.provides.build
.PHONY : blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.provides

blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.provides.build: blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o


# Object files for target zblat1
zblat1_OBJECTS = \
"CMakeFiles/zblat1.dir/zblat1.f.o"

# External object files for target zblat1
zblat1_EXTERNAL_OBJECTS =

blas/testing/zblat1: blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o
blas/testing/zblat1: blas/testing/CMakeFiles/zblat1.dir/build.make
blas/testing/zblat1: blas/libeigen_blas.so
blas/testing/zblat1: blas/testing/CMakeFiles/zblat1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable zblat1"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zblat1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
blas/testing/CMakeFiles/zblat1.dir/build: blas/testing/zblat1

.PHONY : blas/testing/CMakeFiles/zblat1.dir/build

blas/testing/CMakeFiles/zblat1.dir/requires: blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o.requires

.PHONY : blas/testing/CMakeFiles/zblat1.dir/requires

blas/testing/CMakeFiles/zblat1.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing && $(CMAKE_COMMAND) -P CMakeFiles/zblat1.dir/cmake_clean.cmake
.PHONY : blas/testing/CMakeFiles/zblat1.dir/clean

blas/testing/CMakeFiles/zblat1.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/blas/testing /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/blas/testing/CMakeFiles/zblat1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : blas/testing/CMakeFiles/zblat1.dir/depend

