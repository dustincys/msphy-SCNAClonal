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
include test/CMakeFiles/stdlist_4.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/stdlist_4.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/stdlist_4.dir/flags.make

test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o: test/CMakeFiles/stdlist_4.dir/flags.make
test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o: ../test/stdlist.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/stdlist_4.dir/stdlist.cpp.o -c /media/d/github/phy-SCNAClonal/dependencies/eigen/test/stdlist.cpp

test/CMakeFiles/stdlist_4.dir/stdlist.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/stdlist_4.dir/stdlist.cpp.i"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/d/github/phy-SCNAClonal/dependencies/eigen/test/stdlist.cpp > CMakeFiles/stdlist_4.dir/stdlist.cpp.i

test/CMakeFiles/stdlist_4.dir/stdlist.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/stdlist_4.dir/stdlist.cpp.s"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/d/github/phy-SCNAClonal/dependencies/eigen/test/stdlist.cpp -o CMakeFiles/stdlist_4.dir/stdlist.cpp.s

test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.requires:

.PHONY : test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.requires

test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.provides: test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/stdlist_4.dir/build.make test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.provides.build
.PHONY : test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.provides

test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.provides.build: test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o


# Object files for target stdlist_4
stdlist_4_OBJECTS = \
"CMakeFiles/stdlist_4.dir/stdlist.cpp.o"

# External object files for target stdlist_4
stdlist_4_EXTERNAL_OBJECTS =

test/stdlist_4: test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o
test/stdlist_4: test/CMakeFiles/stdlist_4.dir/build.make
test/stdlist_4: test/CMakeFiles/stdlist_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable stdlist_4"
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/stdlist_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/stdlist_4.dir/build: test/stdlist_4

.PHONY : test/CMakeFiles/stdlist_4.dir/build

test/CMakeFiles/stdlist_4.dir/requires: test/CMakeFiles/stdlist_4.dir/stdlist.cpp.o.requires

.PHONY : test/CMakeFiles/stdlist_4.dir/requires

test/CMakeFiles/stdlist_4.dir/clean:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test && $(CMAKE_COMMAND) -P CMakeFiles/stdlist_4.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/stdlist_4.dir/clean

test/CMakeFiles/stdlist_4.dir/depend:
	cd /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/d/github/phy-SCNAClonal/dependencies/eigen /media/d/github/phy-SCNAClonal/dependencies/eigen/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test /media/d/github/phy-SCNAClonal/dependencies/eigen/build_dir/test/CMakeFiles/stdlist_4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/stdlist_4.dir/depend

