%
% ROUTINE DEFINE_PATHS: DEFINES PATHS TO REQUIRED SUBROUTINES
%
% Path to Java classes should be modified according to user's system.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

% turn the warning 'Java class already in static path' off
warning('off','MATLAB:javaclasspath:jarAlreadySpecified');

% PATH TO MATLAB SUBROUTINES
addpath('src')
addpath('src/src_ip4di')
addpath('src/src_nogui_ip4di')
addpath('src/src_mesh2d')
addpath('src/src_IGI')
addpath('cmaps')

% PATH TO JAVA CLASSES FOR IMAGE-GUIDED INVERSION
javaaddpath('src/java-src')

% path to Dave Hale's JTK and IDH libraries
% (to be changed according to user's system)
LIB_PATH='/usr/local/lib/';
javaaddpath([LIB_PATH 'idh/bench/build/classes'])
javaaddpath([LIB_PATH 'jtk/build/libs/edu_mines_jtk.jar'])
javaaddpath([LIB_PATH 'jtk/libs/arpack-java.jar'])
javaaddpath([LIB_PATH 'jtk/libs/netlib-java.jar'])
javaaddpath([LIB_PATH 'jtk/libs/gluegen-rt.jar'])
javaaddpath([LIB_PATH 'jtk/libs/jogl-all.jar'])
javaaddpath([LIB_PATH 'jtk/libs/junit.jar'])
javaaddpath([LIB_PATH 'jtk/libs/jythonlib.jar'])
