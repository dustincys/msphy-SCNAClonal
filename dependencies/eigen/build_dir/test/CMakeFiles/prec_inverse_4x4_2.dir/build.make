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
include test/CMakeFiles/prec_inverse_4x4_2.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/prec_inverse_4x4_2.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/prec_inverse_4x4_2.dir/flags.make

test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o: test/CMakeFiles/prec_inverse_4x4_2.dir/flags.make
test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o: ../test/prec_inverse_4x4.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o -c /media/d/github/phy-SCNAClonal/dependencies/eigen/test/prec_inverse_4x4.cpp

test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.i"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/d/github/phy-SCNAClonal/dependencies/eigen/test/prec_inverse_4x4.cpp > CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.i

test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.s"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/d/github/phy-SCNAClonal/dependencies/eigen/test/prec_inverse_4x4.cpp -o CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.s

test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.requires:

.PHONY : test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.requires

test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.provides: test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/prec_inverse_4x4_2.dir/build.make test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.provides.build
.PHONY : test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.provides

test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.provides.build: test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o


# Object files for target prec_inverse_4x4_2
prec_inverse_4x4_2_OBJECTS = \
"CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o"

# External object files for target prec_inverse_4x4_2
prec_inverse_4x4_2_EXTERNAL_OBJECTS =

test/prec_inverse_4x4_2: test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o
test/prec_inverse_4x4_2: test/CMakeFiles/prec_inverse_4x4_2.dir/build.make
test/prec_inverse_4x4_2: test/CMakeFiles/prec_inverse_4x4_2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable prec_inverse_4x4_2"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/prec_inverse_4x4_2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/prec_inverse_4x4_2.dir/build: test/prec_inverse_4x4_2

.PHONY : test/CMakeFiles/prec_inverse_4x4_2.dir/build

test/CMakeFiles/prec_inverse_4x4_2.dir/requires: test/CMakeFiles/prec_inverse_4x4_2.dir/prec_inverse_4x4.cpp.o.requires

.PHONY : test/CMakeFiles/prec_inverse_4x4_2.dir/requires

test/CMakeFiles/prec_inverse_4x4_2.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -P CMakeFiles/prec_inverse_4x4_2.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/prec_inverse_4x4_2.dir/clean

test/CMakeFiles/prec_inverse_4x4_2.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test/CMakeFiles/prec_inverse_4x4_2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/prec_inverse_4x4_2.dir/depend

