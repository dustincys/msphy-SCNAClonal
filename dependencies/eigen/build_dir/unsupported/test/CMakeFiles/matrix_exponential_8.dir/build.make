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
include unsupported/test/CMakeFiles/matrix_exponential_8.dir/depend.make

# Include the progress variables for this target.
include unsupported/test/CMakeFiles/matrix_exponential_8.dir/progress.make

# Include the compile flags for this target's objects.
include unsupported/test/CMakeFiles/matrix_exponential_8.dir/flags.make

unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o: unsupported/test/CMakeFiles/matrix_exponential_8.dir/flags.make
unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o: ../unsupported/test/matrix_exponential.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o -c /media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/test/matrix_exponential.cpp

unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.i"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/test/matrix_exponential.cpp > CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.i

unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.s"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/test/matrix_exponential.cpp -o CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.s

unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.requires:

.PHONY : unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.requires

unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.provides: unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.requires
	$(MAKE) -f unsupported/test/CMakeFiles/matrix_exponential_8.dir/build.make unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.provides.build
.PHONY : unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.provides

unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.provides.build: unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o


# Object files for target matrix_exponential_8
matrix_exponential_8_OBJECTS = \
"CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o"

# External object files for target matrix_exponential_8
matrix_exponential_8_EXTERNAL_OBJECTS =

unsupported/test/matrix_exponential_8: unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o
unsupported/test/matrix_exponential_8: unsupported/test/CMakeFiles/matrix_exponential_8.dir/build.make
unsupported/test/matrix_exponential_8: unsupported/test/CMakeFiles/matrix_exponential_8.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable matrix_exponential_8"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix_exponential_8.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unsupported/test/CMakeFiles/matrix_exponential_8.dir/build: unsupported/test/matrix_exponential_8

.PHONY : unsupported/test/CMakeFiles/matrix_exponential_8.dir/build

unsupported/test/CMakeFiles/matrix_exponential_8.dir/requires: unsupported/test/CMakeFiles/matrix_exponential_8.dir/matrix_exponential.cpp.o.requires

.PHONY : unsupported/test/CMakeFiles/matrix_exponential_8.dir/requires

unsupported/test/CMakeFiles/matrix_exponential_8.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test && $(CMAKE_COMMAND) -P CMakeFiles/matrix_exponential_8.dir/cmake_clean.cmake
.PHONY : unsupported/test/CMakeFiles/matrix_exponential_8.dir/clean

unsupported/test/CMakeFiles/matrix_exponential_8.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/unsupported/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/unsupported/test/CMakeFiles/matrix_exponential_8.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unsupported/test/CMakeFiles/matrix_exponential_8.dir/depend

