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
include test/CMakeFiles/unalignedassert_3.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/unalignedassert_3.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/unalignedassert_3.dir/flags.make

test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o: test/CMakeFiles/unalignedassert_3.dir/flags.make
test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o: ../test/unalignedassert.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o -c /media/d/github/phy-SCNAClonal/dependencies/eigen/test/unalignedassert.cpp

test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.i"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/d/github/phy-SCNAClonal/dependencies/eigen/test/unalignedassert.cpp > CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.i

test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.s"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/d/github/phy-SCNAClonal/dependencies/eigen/test/unalignedassert.cpp -o CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.s

test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.requires:

.PHONY : test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.requires

test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.provides: test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/unalignedassert_3.dir/build.make test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.provides.build
.PHONY : test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.provides

test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.provides.build: test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o


# Object files for target unalignedassert_3
unalignedassert_3_OBJECTS = \
"CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o"

# External object files for target unalignedassert_3
unalignedassert_3_EXTERNAL_OBJECTS =

test/unalignedassert_3: test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o
test/unalignedassert_3: test/CMakeFiles/unalignedassert_3.dir/build.make
test/unalignedassert_3: test/CMakeFiles/unalignedassert_3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable unalignedassert_3"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/unalignedassert_3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/unalignedassert_3.dir/build: test/unalignedassert_3

.PHONY : test/CMakeFiles/unalignedassert_3.dir/build

test/CMakeFiles/unalignedassert_3.dir/requires: test/CMakeFiles/unalignedassert_3.dir/unalignedassert.cpp.o.requires

.PHONY : test/CMakeFiles/unalignedassert_3.dir/requires

test/CMakeFiles/unalignedassert_3.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -P CMakeFiles/unalignedassert_3.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/unalignedassert_3.dir/clean

test/CMakeFiles/unalignedassert_3.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test/CMakeFiles/unalignedassert_3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/unalignedassert_3.dir/depend

