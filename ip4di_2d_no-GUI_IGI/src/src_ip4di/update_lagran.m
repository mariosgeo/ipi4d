%
% FUNCTION UPDATE_LAGRAN: UPDATE THE LAGRANGIAN MULTIPLIER
%   which controls the weight of the model term
%   in the objective function.
%
% Author: Marios Karaoulis, Colorado School of Mines
%
% Modified by Francois Lavoue', Colorado School of Mines, September 2015.
% Modifications include:
% 1) some explicit comments for more readability,
% 2) more options in the change of the Lagrangian multiplier.
%    The user can now choose either to
%    - keep it constant during inversion
%      (input.decrease_lagrn_flag=0 in 'inversion_parameters.m')
%    - decrease it by a constant factor input.lagrn_reduction
%      (input.decrease_lagrn_flag=1, input.decrease_lagrn_reduction_flag=0)
%    - or decrease it with MK's prescribed factors
%      (input.decrease_lagrn_flag=1, input.decrease_lagrn_reduction_flag=1)

function [input]=update_lagran(itr,ip_cnt,input)
 
%---------------------LAGRANGIAN VALUE AND REDUCTIONS----------------------   

 % Keep original values for ip inversion  
 if (ip_cnt==1 && itr==1) 
 %if not IP
    input.original_lagrn_min=input.lagrn_min; % ACB
    input.original_lagrn_max=input.lagrn_max; %ACB
    input.original_lagrn=input.lagrn;   %Classic
 end

 if (ip_cnt==2 && itr==1) 
 % if IP
    input.lagrn_min=input.original_lagrn_min;
    input.lagrn_max=input.original_lagrn_max;
    input.lagrn=input.original_lagrn; 
 end


 if itr>=1 && input.decrease_lagrn_reduction_flag==1
 % decrease Lagrangian rate of decrease...
 %FL: what the point of hard-coding lagrn_reduction?
    if itr==1 ;input.lagrn_reduction=2;end
    if itr==2 ;input.lagrn_reduction=1.75; end
    if itr==3 ;input.lagrn_reduction=1.5; end
    if itr==4 ;input.lagrn_reduction=1.25; end
    if itr>4  ;input.lagrn_reduction=1; end
 end 


 if itr>=1 && input.decrease_lagrn_flag==1
 % decrease Lagrangian multiplier,
 % ie weight of model term in the objective function
    %%FL: do NOT decrease lagrn_min and lagrn_max in ACB distribution,
    %     just the global Lagrangian weight.
    %input.lagrn_min=input.lagrn_min/input.lagrn_reduction;
    %input.lagrn_max=input.lagrn_max/input.lagrn_reduction;
    input.lagrn=input.lagrn/input.lagrn_reduction;
 end


end   %end function update_lagran


