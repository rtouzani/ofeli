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
include demos/fluid/tins2/CMakeFiles/tins2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include demos/fluid/tins2/CMakeFiles/tins2.dir/compiler_depend.make

# Include the progress variables for this target.
include demos/fluid/tins2/CMakeFiles/tins2.dir/progress.make

# Include the compile flags for this target's objects.
include demos/fluid/tins2/CMakeFiles/tins2.dir/flags.make

demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.o: demos/fluid/tins2/CMakeFiles/tins2.dir/flags.make
demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.o: ../demos/fluid/tins2/main.cpp
demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.o: demos/fluid/tins2/CMakeFiles/tins2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/touzani/ofeli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.o"
	cd /Users/touzani/ofeli/build/demos/fluid/tins2 && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.o -MF CMakeFiles/tins2.dir/main.cpp.o.d -o CMakeFiles/tins2.dir/main.cpp.o -c /Users/touzani/ofeli/demos/fluid/tins2/main.cpp

demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tins2.dir/main.cpp.i"
	cd /Users/touzani/ofeli/build/demos/fluid/tins2 && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/touzani/ofeli/demos/fluid/tins2/main.cpp > CMakeFiles/tins2.dir/main.cpp.i

demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tins2.dir/main.cpp.s"
	cd /Users/touzani/ofeli/build/demos/fluid/tins2 && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/touzani/ofeli/demos/fluid/tins2/main.cpp -o CMakeFiles/tins2.dir/main.cpp.s

# Object files for target tins2
tins2_OBJECTS = \
"CMakeFiles/tins2.dir/main.cpp.o"

# External object files for target tins2
tins2_EXTERNAL_OBJECTS =

demos/fluid/tins2/tins2: demos/fluid/tins2/CMakeFiles/tins2.dir/main.cpp.o
demos/fluid/tins2/tins2: demos/fluid/tins2/CMakeFiles/tins2.dir/build.make
demos/fluid/tins2/tins2: src/libofeli.a
demos/fluid/tins2/tins2: demos/fluid/tins2/CMakeFiles/tins2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/touzani/ofeli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tins2"
	cd /Users/touzani/ofeli/build/demos/fluid/tins2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tins2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/fluid/tins2/CMakeFiles/tins2.dir/build: demos/fluid/tins2/tins2
.PHONY : demos/fluid/tins2/CMakeFiles/tins2.dir/build

demos/fluid/tins2/CMakeFiles/tins2.dir/clean:
	cd /Users/touzani/ofeli/build/demos/fluid/tins2 && $(CMAKE_COMMAND) -P CMakeFiles/tins2.dir/cmake_clean.cmake
.PHONY : demos/fluid/tins2/CMakeFiles/tins2.dir/clean

demos/fluid/tins2/CMakeFiles/tins2.dir/depend:
	cd /Users/touzani/ofeli/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/touzani/ofeli /Users/touzani/ofeli/demos/fluid/tins2 /Users/touzani/ofeli/build /Users/touzani/ofeli/build/demos/fluid/tins2 /Users/touzani/ofeli/build/demos/fluid/tins2/CMakeFiles/tins2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/fluid/tins2/CMakeFiles/tins2.dir/depend

