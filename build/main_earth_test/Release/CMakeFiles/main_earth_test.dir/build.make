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
CMAKE_BINARY_DIR = /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release

# Include any dependencies generated for this target.
include CMakeFiles/main_earth_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main_earth_test.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main_earth_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main_earth_test.dir/flags.make

CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o: CMakeFiles/main_earth_test.dir/flags.make
CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o: ../../../main_earth_test.cpp
CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o: CMakeFiles/main_earth_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o -MF CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o.d -o CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o -c /Users/ben/Desktop/allo/allolib_playground/ocean/main_earth_test.cpp

CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ben/Desktop/allo/allolib_playground/ocean/main_earth_test.cpp > CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.i

CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ben/Desktop/allo/allolib_playground/ocean/main_earth_test.cpp -o CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.s

# Object files for target main_earth_test
main_earth_test_OBJECTS = \
"CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o"

# External object files for target main_earth_test
main_earth_test_EXTERNAL_OBJECTS =

../../../bin/main_earth_test: CMakeFiles/main_earth_test.dir/ocean/main_earth_test.cpp.o
../../../bin/main_earth_test: CMakeFiles/main_earth_test.dir/build.make
../../../bin/main_earth_test: ../../../../allolib/build/Release/libal.a
../../../bin/main_earth_test: ../../../../allolib/build/al_ext/assets3d/libal_assets3d.a
../../../bin/main_earth_test: ../../../../allolib/build/al_ext/openvr/libal_openvr.a
../../../bin/main_earth_test: ../../../../allolib/build/al_ext/soundfile/libal_soundfile.a
../../../bin/main_earth_test: ../../../../allolib/build/al_ext/statedistribution/libal_statedistribution.a
../../../bin/main_earth_test: /opt/homebrew/lib/libassimp.dylib
../../../bin/main_earth_test: ../../../../allolib/build/Release/libal.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/rtaudio/librtaudio.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/Gamma/lib/libGamma.a
../../../bin/main_earth_test: /opt/homebrew/lib/libsndfile.dylib
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/glfw/src/libglfw3.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/glad/libglad.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/rtmidi/librtmidi.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/libimgui.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/liboscpack.a
../../../bin/main_earth_test: ../../../../allolib/build/Release/external/libserial.a
../../../bin/main_earth_test: CMakeFiles/main_earth_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../bin/main_earth_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main_earth_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main_earth_test.dir/build: ../../../bin/main_earth_test
.PHONY : CMakeFiles/main_earth_test.dir/build

CMakeFiles/main_earth_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main_earth_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main_earth_test.dir/clean

CMakeFiles/main_earth_test.dir/depend:
	cd /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ben/Desktop/allo/allolib_playground /Users/ben/Desktop/allo/allolib_playground /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release /Users/ben/Desktop/allo/allolib_playground/ocean/build/main_earth_test/Release/CMakeFiles/main_earth_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main_earth_test.dir/depend

