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

# Utility rule file for demos.

# Include the progress variables for this target.
include demos/CMakeFiles/demos.dir/progress.make

demos: demos/CMakeFiles/demos.dir/build.make

.PHONY : demos

# Rule to build all files generated by this target.
demos/CMakeFiles/demos.dir/build: demos

.PHONY : demos/CMakeFiles/demos.dir/build

demos/CMakeFiles/demos.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/demos && $(CMAKE_COMMAND) -P CMakeFiles/demos.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/demos.dir/clean

demos/CMakeFiles/demos.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/demos /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/demos /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/demos/CMakeFiles/demos.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/demos.dir/depend

