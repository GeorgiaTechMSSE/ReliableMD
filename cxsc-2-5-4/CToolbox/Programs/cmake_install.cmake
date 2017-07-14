# Install script for directory: /home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/xev_ex2.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/rpe_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/hess_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/jac_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/fastlss_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/gop1_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/lss_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/xev_ex3.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/ddf_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/xev_ex1.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/gop_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/cpz_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/lop_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/nlss_ex.cpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/nlfz_ex.cpp"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/hess_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/hess_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/cpz_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/cpz_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/ddf_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/ddf_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/fastlss_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fastlss_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/gop1_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop1_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/gop_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/gop_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/jac_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/jac_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/lop_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lop_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/lss_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/lss_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/nlfz_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlfz_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/nlss_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/nlss_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/rpe_ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/rpe_ex")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/xev_ex1")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex1")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/xev_ex2")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex2")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3"
         RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/xev_ex3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3"
         OLD_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib:/home/anhvt89/Documents/c++lib/cxsc-2-5-4:"
         NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/xev_ex3")
    endif()
  endif()
endif()

