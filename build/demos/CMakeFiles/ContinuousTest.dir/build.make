# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/touzani/ofeli

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/touzani/ofeli/build

# Utility rule file for ContinuousTest.

# Include any custom commands dependencies for this target.
include demos/CMakeFiles/ContinuousTest.dir/compiler_depend.make

# Include the progress variables for this target.
include demos/CMakeFiles/ContinuousTest.dir/progress.make

demos/CMakeFiles/ContinuousTest:
	cd /Users/touzani/ofeli/build/demos && /opt/local/bin/ctest -D ContinuousTest

ContinuousTest: demos/CMakeFiles/ContinuousTest
ContinuousTest: demos/CMakeFiles/ContinuousTest.dir/build.make
.PHONY : ContinuousTest

# Rule to build all files generated by this target.
demos/CMakeFiles/ContinuousTest.dir/build: ContinuousTest
.PHONY : demos/CMakeFiles/ContinuousTest.dir/build

demos/CMakeFiles/ContinuousTest.dir/clean:
	cd /Users/touzani/ofeli/build/demos && $(CMAKE_COMMAND) -P CMakeFiles/ContinuousTest.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/ContinuousTest.dir/clean

demos/CMakeFiles/ContinuousTest.dir/depend:
	cd /Users/touzani/ofeli/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/touzani/ofeli /Users/touzani/ofeli/demos /Users/touzani/ofeli/build /Users/touzani/ofeli/build/demos /Users/touzani/ofeli/build/demos/CMakeFiles/ContinuousTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/ContinuousTest.dir/depend

