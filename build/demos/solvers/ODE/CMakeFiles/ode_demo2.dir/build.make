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
include demos/solvers/ODE/CMakeFiles/ode_demo2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include demos/solvers/ODE/CMakeFiles/ode_demo2.dir/compiler_depend.make

# Include the progress variables for this target.
include demos/solvers/ODE/CMakeFiles/ode_demo2.dir/progress.make

# Include the compile flags for this target's objects.
include demos/solvers/ODE/CMakeFiles/ode_demo2.dir/flags.make

demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o: demos/solvers/ODE/CMakeFiles/ode_demo2.dir/flags.make
demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o: ../demos/solvers/ODE/ode_demo2.cpp
demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o: demos/solvers/ODE/CMakeFiles/ode_demo2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/touzani/ofeli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o"
	cd /Users/touzani/ofeli/build/demos/solvers/ODE && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o -MF CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o.d -o CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o -c /Users/touzani/ofeli/demos/solvers/ODE/ode_demo2.cpp

demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ode_demo2.dir/ode_demo2.cpp.i"
	cd /Users/touzani/ofeli/build/demos/solvers/ODE && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/touzani/ofeli/demos/solvers/ODE/ode_demo2.cpp > CMakeFiles/ode_demo2.dir/ode_demo2.cpp.i

demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ode_demo2.dir/ode_demo2.cpp.s"
	cd /Users/touzani/ofeli/build/demos/solvers/ODE && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/touzani/ofeli/demos/solvers/ODE/ode_demo2.cpp -o CMakeFiles/ode_demo2.dir/ode_demo2.cpp.s

# Object files for target ode_demo2
ode_demo2_OBJECTS = \
"CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o"

# External object files for target ode_demo2
ode_demo2_EXTERNAL_OBJECTS =

demos/solvers/ODE/ode_demo2: demos/solvers/ODE/CMakeFiles/ode_demo2.dir/ode_demo2.cpp.o
demos/solvers/ODE/ode_demo2: demos/solvers/ODE/CMakeFiles/ode_demo2.dir/build.make
demos/solvers/ODE/ode_demo2: src/libofeli.a
demos/solvers/ODE/ode_demo2: demos/solvers/ODE/CMakeFiles/ode_demo2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/touzani/ofeli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ode_demo2"
	cd /Users/touzani/ofeli/build/demos/solvers/ODE && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ode_demo2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/solvers/ODE/CMakeFiles/ode_demo2.dir/build: demos/solvers/ODE/ode_demo2
.PHONY : demos/solvers/ODE/CMakeFiles/ode_demo2.dir/build

demos/solvers/ODE/CMakeFiles/ode_demo2.dir/clean:
	cd /Users/touzani/ofeli/build/demos/solvers/ODE && $(CMAKE_COMMAND) -P CMakeFiles/ode_demo2.dir/cmake_clean.cmake
.PHONY : demos/solvers/ODE/CMakeFiles/ode_demo2.dir/clean

demos/solvers/ODE/CMakeFiles/ode_demo2.dir/depend:
	cd /Users/touzani/ofeli/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/touzani/ofeli /Users/touzani/ofeli/demos/solvers/ODE /Users/touzani/ofeli/build /Users/touzani/ofeli/build/demos/solvers/ODE /Users/touzani/ofeli/build/demos/solvers/ODE/CMakeFiles/ode_demo2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/solvers/ODE/CMakeFiles/ode_demo2.dir/depend

