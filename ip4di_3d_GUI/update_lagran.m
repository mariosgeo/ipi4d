function [input]=update_lagran(itr,ip_cnt,tmp_flag,input)
 
 
 
 
   %---------------------LAGRANGIAN VALUE AND REDUCTIONS----------------------   
 % Keep original values for ip inversion  
if (ip_cnt==1 && itr==1) 
    input.original_lagrn_min=input.lagrn_min; % ACB
    input.original_lagrn_max=input.lagrn_max; %ACB
    input.original_lagrn=input.lagrn;   %Classic
end
if (ip_cnt==2 && itr==1) 
    input.lagrn_min=input.original_lagrn_min;
    input.lagrn_max=input.original_lagrn_max;
    input.lagrn=input.original_lagrn; 
end

if itr>1 && itr<=5 && tmp_flag==1

    if itr==1 ;input.lagrn_reduction=2;end
    if itr==2 ;input.lagrn_reduction=1.75; end
    if itr==3 ;input.lagrn_reduction=1.5; end
    if itr==4 ;input.lagrn_reduction=1.25; end
    if itr>4  ;input.lagrn_reduction=1; end
    
    input.lagrn_min=input.lagrn_min/input.lagrn_reduction;
    input.lagrn_max=input.lagrn_max/input.lagrn_reduction;
    input.lagrn=input.lagrn/input.lagrn_reduction;
end 

end