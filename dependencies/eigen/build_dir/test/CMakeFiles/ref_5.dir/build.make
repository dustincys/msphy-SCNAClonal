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
include test/CMakeFiles/ref_5.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/ref_5.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/ref_5.dir/flags.make

test/CMakeFiles/ref_5.dir/ref.cpp.o: test/CMakeFiles/ref_5.dir/flags.make
test/CMakeFiles/ref_5.dir/ref.cpp.o: ../test/ref.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/ref_5.dir/ref.cpp.o"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ref_5.dir/ref.cpp.o -c /media/d/github/phy-SCNAClonal/dependencies/eigen/test/ref.cpp

test/CMakeFiles/ref_5.dir/ref.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ref_5.dir/ref.cpp.i"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/d/github/phy-SCNAClonal/dependencies/eigen/test/ref.cpp > CMakeFiles/ref_5.dir/ref.cpp.i

test/CMakeFiles/ref_5.dir/ref.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ref_5.dir/ref.cpp.s"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/d/github/phy-SCNAClonal/dependencies/eigen/test/ref.cpp -o CMakeFiles/ref_5.dir/ref.cpp.s

test/CMakeFiles/ref_5.dir/ref.cpp.o.requires:

.PHONY : test/CMakeFiles/ref_5.dir/ref.cpp.o.requires

test/CMakeFiles/ref_5.dir/ref.cpp.o.provides: test/CMakeFiles/ref_5.dir/ref.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/ref_5.dir/build.make test/CMakeFiles/ref_5.dir/ref.cpp.o.provides.build
.PHONY : test/CMakeFiles/ref_5.dir/ref.cpp.o.provides

test/CMakeFiles/ref_5.dir/ref.cpp.o.provides.build: test/CMakeFiles/ref_5.dir/ref.cpp.o


# Object files for target ref_5
ref_5_OBJECTS = \
"CMakeFiles/ref_5.dir/ref.cpp.o"

# External object files for target ref_5
ref_5_EXTERNAL_OBJECTS =

test/ref_5: test/CMakeFiles/ref_5.dir/ref.cpp.o
test/ref_5: test/CMakeFiles/ref_5.dir/build.make
test/ref_5: test/CMakeFiles/ref_5.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ref_5"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ref_5.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/ref_5.dir/build: test/ref_5

.PHONY : test/CMakeFiles/ref_5.dir/build

test/CMakeFiles/ref_5.dir/requires: test/CMakeFiles/ref_5.dir/ref.cpp.o.requires

.PHONY : test/CMakeFiles/ref_5.dir/requires

test/CMakeFiles/ref_5.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -P CMakeFiles/ref_5.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/ref_5.dir/clean

test/CMakeFiles/ref_5.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test/CMakeFiles/ref_5.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/ref_5.dir/depend

