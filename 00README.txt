
 IMAGE-GUIDED IP4DI without Graphical User Interface
    IGI   -   IP4DI   -   no GUI

This is a modified version of M. Karaoulis' IP4DI program (MK, Colorado School of Mines),
derived from Jieyi Zhou's version as of July 23, 2015 (JZ, Colorado School of Mines),
and modified by F. Lavoue' (FL, Colorado School of Mines, francois.lavoue@ens-lyon.org).

FL's modifications mainly concerned the Image-Guided Inversion (IGI) approach, 
using JZ's work and Dave Hale's (DH) Java routines for image processing.
But DH's Java routines require JDK 1.7, which turned out to be incompatible
with MK's IP4DI Graphical User Guide (GUI), that works best with JDK 1.6.
Therefore, significant modifications in the organisation of the code have also
been made in order to use IP4DI routines without the GUI. Besides, not using
the GUI also enables to use the code in a more advanced way ("nohup"-like
execution from a Unix terminal, automatic sensitivity analysis and L-curves
computation, for instance).

Programs also have been heavily commented to facilitate readability and use.
As a consequence, input files and main programs should now be self-explanatory. 
If not, feel free to address your questions to francois.lavoue@ens-lyon.org.

Francois Lavoue', December 8, 2015


----------------------------------------------------------------------------
ORGANISATION OF THE CODE:

This readme file is in the main directory from where the code should be launched.
It is recommended to keep this directory as clean as possible (avoid having too
much files in this main directory, rather organize your files in folders 'data',
'results', 'tools', etc).

Source codes are in ./src, classified in
- ./src/src_ip4di: M. Karaoulis' IP4DI routines (eventually modified by JZ and FL).
- ./src/src_nogui_ip4di: routines derived from IP4DI, for a use without GUI (FL).
- ./src/src_mesh2d: routines for mesh generation (by D. Engwirda, Univ. Sydney, Sept. 2005).
- ./src/src_IGI: J. Zhou's routines for Image-Guided Inversion (modified by FL).
Old source files have been kept in 'old*' directories for development purposes
(track changes and backup).

  The main programs, 'EXEC_FORWARD.m' and 'EXEC_INVERSION.m', have been kept in main
directory, as well as:
- 'define_paths.m' that defines the path to Matlab subroutines and Java classes,
- the input routines 'forward_parameters.m', 'inversion_parameters.m' and 'plot_parameters.m',
- and the 'meta' routine 'sensitivity_analysis.m' which enables to repeat the inversion
  for several values of a given parameter (enabling for example to compute L-curves).

  These input routines are the only files that must be modified to use the code, which is
why they are in the main directory for easy editing (while the files 'EXEC*.m' and 
'define_paths.m' must be present in the main directory to launch the code).

  A number of simple programs have been put in 'tools' directory for plotting
inversion results once the program terminated. These programs are intended to be
launched in the main directory.


----------------------------------------------------------------------------
USING THE CODE

To run a forward simulation, edit the file 'forward_parameters.m', 
launch Matlab from the current directory and type 'EXEC_FORWARD' in the Matlab console.

To perform an inversion, edit the files 'forward_parameters.m' and 'inverse_parameters.m',
launch Matlab from the current directory and type 'EXEC_INVERSION' in the Matlab console.

To perform a sensitivity analysis, edit the files 'forward_parameters.m', 'inverse_parameters.m'
and 'sensitivity_analysis.m' to choose the varying parameter and its range of values,
launch Matlab from the current directory and type 'EXEC_INVERSION' in the Matlab console.

To run the code in a nohup-like way, instead of opening Matlab and launching the code from
Matlab console, type in a Terminal window:
------
matlab -nosplash -nodisplay -nodesktop < EXEC_INVERSION.m >log.out &
------
This enables the user to log out without interrupting the job.


----------------------------------------------------------------------------
USING THE CODE FOR IMAGE-GUIDED INVERSION (IGI): 

To use IGI, your Java version needs to be 1.7 (you must have installed Java SE JDK 7),
and you must set up the following environment variables (in command line or in your .bashrc):
------
export JAVA_HOME=/usr/lib/jvm/jdk1.7.0_40
------
(just an example, this needs to be the folder that contains /bin/java)
and
------
export MATLAB_JAVA=$JAVA_HOME/jre
------
You can check the version of Java used by Matlab by typing in the Matlab console
-------
version -java
-------

Before using IGI, you need to download Mines Java Toolkit (JTK) from
https://github.com/dhale/jtk
(it is more convenient to use git from the terminal to do it.)
Follow the instructions on this page to install it.

Then, also download Dave Hale's package IDH in the same directory as JTK, from 
https://github.com/dhale/idh
and copy the folder 'gradle' and the file 'gradlew' from folder 'jtk' to 'idh'.
Install IDH the same way as JTK.

The directory where you installed JTK and IDH (e.g. /usr/local/lib) must be specified 
as LIB_PATH in 'define_paths.m' for Matlab to find the installed Java classes.

NB: JDK 7 is not only incompatible with original IP4DI's GUI, it might also cause
problems with Matlab GUI and with figure display. Therefore it is recommended to
launch Matlab without using the GUI, directly in a Terminal window by typing
------
matlab -nodesktop
------
and to avoid plotting intermediate results during the inversion (unless for debug/check), 
because plots may slow down the execution of the program (see plot options in
'forward_parameters.m' and 'inversion_parameters.m'). 
Figure can be displayed after inversion with the 'tools' programs, and preferably
in a new Terminal session where the environment variables $JAVA_HOME and $MATLAB_JAVA
are reset to enable the use of JDK 1.6:
------
export JAVA_HOME=
export MATLAB_JAVA=
------


----------------------------------------------------------------------------
INPUT/OUTPUT FORMATS:

- The input models must be specified in an ASCII file of 3 or 4 columns for
        x  z  res.                      % for purely real resistivity
  or    x  z  Re(res.)   Im(res.)       % if input.cmplx_flag=1 is selected
  or    x  z  amp(res.)  phase(res.)    % if input.cmplx_flag=2,
  with as many lines as model parameters.

  NB: if input.interp_model_flag=0, input model must match the mesh, with x and
  z corresponding to the center of the cells. If the mesh is not known (e.g. not
  not built yet), the user can select input.interp_model_flag=1 in forward_parameters.m
  to interpolate an arbitrary input model on the generated mesh. This will output
  the interpolated model in file input.file_model_out for future use.

- Output models are delivered in the same format. If input.print_intermediate_results=1
  is selected in inversion_parameters.m, the models of intermediate iterations are
  appended in additional columns.

- Input data must be specified in an ASCII file of 9 to 10 columns that specify,
  for each data point, the electrode locations and the data value:
       xA zA xB zB xM zM xN zN data_DC                 % DC data
  or   xA zA xB zB xM zM xN zN Re(data)  Im(data)      % SIP data (input.cmplx_flag=1)
  or   xA zA xB zB xM zM xN zN amp(data) phase(data)   % SIP data (input.cmplx_flag=1)
    Values can be specified in terms of resistivities (input.data_type=1) or
  resistances (input.data_type=2).
    For pole-pole or pole-dipole arrays, columns for B and N electrodes must be present
  but their coordinates must be identical for all data points (just put 0).
    IP4DI can also read Res2Dinv data files (set input.res2d_flag=1 in forward_parameters.m).

- Output data are delivered in the same format. If input.print_intermediate_results=1
  is selected in inversion_parameters.m, synthetic data computed at intermediate
  iterations are appended in additional columns.

- After mesh generation, a mesh file is delivered containing the Matlab structure
  of all mesh variables. This file can be read using 'importdata' in Matlab or 
  'load' in Octave.

- If input.output_variables=1 is selected in inversion_parameters.m, inversion also
  delivers Matlab structures 'input', 'mesh' and 'final' at the end of the program.
  The 'final' structure contains inversion results.

- Finally, the input model that serves as a training image for the image-guided
  approach must be provided in a binary file as a matrix of size [nz,nx], with 
  nz the nb of rows (input.n1_TI in inversion_parameters.m) and nx the nb of columns
  (input.n2_TI). The matrix must be ordered column by column, according to Matlab
  convention. The corresponding image is assumed to span in the range 
  [ 0:(nz-1)*h , 0:(nx-1)*h ] in the referential of IP4DI model, with h the grid
  step of the training image (which is generally more finely discretized than the
  reconstructed model to capture the structures of interest).


----------------------------------------------------------------------------
TIPS AND TRICKS

There are a number of tricks to know in order to run the program smoothly. Here are
the most important ones:

- Plots: as said earlier, IGI needs JDK1.7 in order to use DH's JTK library, and this
         makes the display of plots unstable. A solution is to first run IGI with JDK1.7
         (ie with the suitable JAVA_HOME and MATLAB_JAVA in the environment) and to then
         plot the inversion results in a new session after resetting the environment
         variables (see section IGI).

- Mesh: it is recommended to plot the mesh (using input.plot_mesh=1) in order to check
        its consistency. In particular, when using small models (~1m), mesh generation
        may result in weird and useless refinement around some electrodes.


----------------------------------------------------------------------------
MODIFYING THE CODE:

You are free to modify the code and make your own developments, as long as you
respect the legal statements (see 0legal_statement.txt). However, if your
developments aim at being integrated in the stable version of the code on the
long term, or simply at being distributed, some guiding rules should be followed
to ensure a constructive collaborative work, and avoid ending up with too many
messy code versions:
- if possible, consult the current authors before and during your own developments,
  to facilitate their long-term integration.
- add your developments in separate subroutines and call them within the pre-existing 
  workflow. This make the identification of the (new) tasks easier.
- (extensively) comment your code to make it easily understandable. Eventually 
  describe the new functionalities in this readme file.
- give to your variables understandable names (rather 'grad', 'cov' than 'tmp1', 'tmp2').


----------------------------------------------------------------------------
CURRENT BUGS AND FUTURE MODIFICATIONS:

- /!\ Time-lapse inversion still in progress:
      - parameter update in kim_inversion2.m should be validated (a model term is missing).
      - compatibility with ACB (Active Constrained Balancing) must be implemented.
- /!\ ACT (Active Time Constraints) not validated
      (the GUI routine preprocessing.m should be adapted for use without the GUI).
- /!\ IP inversion not validated (should be OK but there might be bugs concerning
      the inversion of the 2nd component, ip_cnt=2).

- make compatible with Octave:
  - avoid reference to importdata function (use 'load' instead),
  - find an alternative to the obsolete tsearch function
    for the generation of the mesh (mesh2d/meshfaces.m, l.199).


----------------------------------------------------------------------------
LOG OF PAST MODIFICATIONS:

- v1: March 2, 2016 (most recent version)

  - correct src/src_ip4di/info_save.m:
    - remove useless loop over nb of param and data
    in the printing of intermediate results.
    - correct output format of final.res_param1_vs_it (real(it1) imag(it1) real(it2) imag(it2),
                                            and not    real(it1) real(it2) imag(it1) imag(it2)).

  - introduce post-inversion image-guided interpolation of the final model
    (both as the subroutine 'src/src_IGI/image_guided_interpolation.m'
     and as the self-standing program 'tools/post_inversion_image_guided_interpolation.m').

- v0: October 25, 2015.

