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
CMAKE_SOURCE_DIR = /home/anhvt89/Documents/c++lib/cxsc-2-5-4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/anhvt89/Documents/c++lib/cxsc-2-5-4

# Include any dependencies generated for this target.
include CToolbox/test/test_programms/CMakeFiles/lss_test.dir/depend.make

# Include the progress variables for this target.
include CToolbox/test/test_programms/CMakeFiles/lss_test.dir/progress.make

# Include the compile flags for this target's objects.
include CToolbox/test/test_programms/CMakeFiles/lss_test.dir/flags.make

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/flags.make
CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o: CToolbox/test/test_programms/lss_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o"
	cd /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lss_test.dir/lss_test.o -c /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms/lss_test.cpp

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lss_test.dir/lss_test.i"
	cd /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms/lss_test.cpp > CMakeFiles/lss_test.dir/lss_test.i

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lss_test.dir/lss_test.s"
	cd /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms/lss_test.cpp -o CMakeFiles/lss_test.dir/lss_test.s

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.requires:

.PHONY : CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.requires

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.provides: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.requires
	$(MAKE) -f CToolbox/test/test_programms/CMakeFiles/lss_test.dir/build.make CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.provides.build
.PHONY : CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.provides

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.provides.build: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o


# Object files for target lss_test
lss_test_OBJECTS = \
"CMakeFiles/lss_test.dir/lss_test.o"

# External object files for target lss_test
lss_test_EXTERNAL_OBJECTS =

CToolbox/test/test_programms/lss_test: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o
CToolbox/test/test_programms/lss_test: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/build.make
CToolbox/test/test_programms/lss_test: libcxsc.so.2.5.4
CToolbox/test/test_programms/lss_test: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lss_test"
	cd /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lss_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CToolbox/test/test_programms/CMakeFiles/lss_test.dir/build: CToolbox/test/test_programms/lss_test

.PHONY : CToolbox/test/test_programms/CMakeFiles/lss_test.dir/build

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/requires: CToolbox/test/test_programms/CMakeFiles/lss_test.dir/lss_test.o.requires

.PHONY : CToolbox/test/test_programms/CMakeFiles/lss_test.dir/requires

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/clean:
	cd /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms && $(CMAKE_COMMAND) -P CMakeFiles/lss_test.dir/cmake_clean.cmake
.PHONY : CToolbox/test/test_programms/CMakeFiles/lss_test.dir/clean

CToolbox/test/test_programms/CMakeFiles/lss_test.dir/depend:
	cd /home/anhvt89/Documents/c++lib/cxsc-2-5-4 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/anhvt89/Documents/c++lib/cxsc-2-5-4 /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms /home/anhvt89/Documents/c++lib/cxsc-2-5-4 /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms/CMakeFiles/lss_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CToolbox/test/test_programms/CMakeFiles/lss_test.dir/depend

