# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.21.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.21.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ben/Desktop/allo/allolib_playground

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth/Release

# Utility rule file for main_earth_run_debug.

# Include any custom commands dependencies for this target.
include CMakeFiles/main_earth_run_debug.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main_earth_run_debug.dir/progress.make

CMakeFiles/main_earth_run_debug: ../../../bin/main_earth
	cd /Users/ben/Desktop/allo/allolib_playground/ocean/bin && lldb -o\ run ./main_earthd

main_earth_run_debug: CMakeFiles/main_earth_run_debug
main_earth_run_debug: CMakeFiles/main_earth_run_debug.dir/build.make
.PHONY : main_earth_run_debug

# Rule to build all files generated by this target.
CMakeFiles/main_earth_run_debug.dir/build: main_earth_run_debug
.PHONY : CMakeFiles/main_earth_run_debug.dir/build

CMakeFiles/main_earth_run_debug.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main_earth_run_debug.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main_earth_run_debug.dir/clean

CMakeFiles/main_earth_run_debug.dir/depend:
	cd /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth/Release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ben/Desktop/allo/allolib_playground /Users/ben/Desktop/allo/allolib_playground /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth/Release /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth/Release /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth/Release/CMakeFiles/main_earth_run_debug.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main_earth_run_debug.dir/depend

