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

# Include any dependencies generated for this target.
include demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/compiler_depend.make

# Include the progress variables for this target.
include demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/progress.make

# Include the compile flags for this target's objects.
include demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/flags.make

demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o: demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/flags.make
demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o: ../demos/solvers/NLAS/nl_demo1.cpp
demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o: demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/touzani/ofeli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o"
	cd /Users/touzani/ofeli/build/demos/solvers/NLAS && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o -MF CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o.d -o CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o -c /Users/touzani/ofeli/demos/solvers/NLAS/nl_demo1.cpp

demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nl_demo1.dir/nl_demo1.cpp.i"
	cd /Users/touzani/ofeli/build/demos/solvers/NLAS && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/touzani/ofeli/demos/solvers/NLAS/nl_demo1.cpp > CMakeFiles/nl_demo1.dir/nl_demo1.cpp.i

demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nl_demo1.dir/nl_demo1.cpp.s"
	cd /Users/touzani/ofeli/build/demos/solvers/NLAS && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/touzani/ofeli/demos/solvers/NLAS/nl_demo1.cpp -o CMakeFiles/nl_demo1.dir/nl_demo1.cpp.s

# Object files for target nl_demo1
nl_demo1_OBJECTS = \
"CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o"

# External object files for target nl_demo1
nl_demo1_EXTERNAL_OBJECTS =

demos/solvers/NLAS/nl_demo1: demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/nl_demo1.cpp.o
demos/solvers/NLAS/nl_demo1: demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/build.make
demos/solvers/NLAS/nl_demo1: src/libofeli.a
demos/solvers/NLAS/nl_demo1: demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/touzani/ofeli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable nl_demo1"
	cd /Users/touzani/ofeli/build/demos/solvers/NLAS && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nl_demo1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/build: demos/solvers/NLAS/nl_demo1
.PHONY : demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/build

demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/clean:
	cd /Users/touzani/ofeli/build/demos/solvers/NLAS && $(CMAKE_COMMAND) -P CMakeFiles/nl_demo1.dir/cmake_clean.cmake
.PHONY : demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/clean

demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/depend:
	cd /Users/touzani/ofeli/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/touzani/ofeli /Users/touzani/ofeli/demos/solvers/NLAS /Users/touzani/ofeli/build /Users/touzani/ofeli/build/demos/solvers/NLAS /Users/touzani/ofeli/build/demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/solvers/NLAS/CMakeFiles/nl_demo1.dir/depend

