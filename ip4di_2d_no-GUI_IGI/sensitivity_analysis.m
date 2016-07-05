%
% Function SENSITIVITY_ANALYSIS: repeat inversion for a range of input parameter values.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function [input,mesh,fem,final]=sensitivity_analysis(input,mesh)

 disp(' ')
 disp('===================================')
 disp('=   PERFORM SENSITIVITY ANALYSIS  =')
 disp('===================================')
 disp(' ') 

 % define range for varying parameter
 % (here Lagrangian multiplier for drawing a L-curve)
 parameter_name='lagrn'
 parameter_range=[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1];

 %
 % LOOP OVER PARAMETER
 %
 for ipar=1:length(parameter_range)

     par_i=parameter_range(ipar);

     % /!\ this line should be changed
     %     according to chosen parameter
     input.lagrn=par_i;

     %
     % OUTPUT FILE NAMES
     % (as many outputs as tested parameter values)
     %

     % file name for info vs. it
     input.file_info=['info_vs_it_' parameter_name num2str(par_i) '.txt'];

     % file names for final data and models
     input.file_data_inv_out=['data_final_' parameter_name num2str(par_i) '.dat'];
     input.file_model_out=['model_final_' parameter_name num2str(par_i) '.dat'];
     input.file_model_interp=['model_final_interp_' parameter_name num2str(par_i) '.dat'];

     % file names for intermediate data and models
     input.print_intermediate_results=1;   % 1->yes, 0->no
     input.file_data_inv_inter=['data_inter_' parameter_name num2str(par_i) '.dat'];
     input.file_model_inter=['model_inter_' parameter_name num2str(par_i) '.dat'];

     % file names for final Matlab structures,
     % enabling to find any variable again
     input.output_variables=1;   % 1->save structures, 0->don't to avoid storing large files
     input.file_input_struct=['struct_input_' parameter_name num2str(par_i) '.mat'];
     input.file_mesh_struct=['struct_mesh_' parameter_name num2str(par_i) '.mat'];
     %input.file_fem_struct=['struct_fem_' parameter_name num2str(par_i) '.mat'];   %usually big and not very useful
     input.file_final_struct=['struct_final_' parameter_name num2str(par_i) '.mat'];

     % RUN INVERSION
     [fem,final]=main(input,mesh);

 end   %end loop over parameter of interest

end   %end function sensitivity_analysis
