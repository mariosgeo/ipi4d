function main(input,mesh)
ex=0;

if matlabpool('size') == 0
       matlabpool open
end
        
global textt
fem=[]; % initialize fem as empty matrix
final=[];

mesh=smooth_mtx_surface4(input,mesh);set(textt,'String','********* SMOOTH MATRIX ********');drawnow update;


% Calculate k and g
[mesh.k,mesh.g]=rwaven(mesh.probe_spacing,mesh.probe_spacing*mesh.max_n);
% rdis=double([mesh.probe_spacing:mesh.probe_spacing:mesh.max_n*mesh.probe_spacing]);
% k_init= logspace(-2,0.3,5);
% [mesh.k,mesh.g,obj,err]=xu_inversion(rdis,k_init);

set(textt,'String','****** INVERSION  STARTS *****'); drawnow update;

if input.time_lapse_flag==0
    for ip_cnt=1:input.ip_num
        for itr=1:input.itn+1
            fem=mes_control_fast(itr,input,mesh,fem,0);
            [ex,fem,mesh]=rms(itr,ip_cnt,input,mesh,fem);
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);
            if (ex==1 || ex==2) ;final.itr=itr-1;break; end
            [input]=update_lagran(itr,ip_cnt,1,input);
            [mesh,fem,input]=invert_cntr(itr,1,ip_cnt,input,mesh,fem);
            [fem,mesh,input]=prop_change(itr,input,mesh,fem);
            auto_contour(1,itr,ip_cnt,input,mesh,fem);
        end
        if itr==input.itn+1; final.itr=itr;end
        if input.ip_flag==1
            [input,mesh]=ip_calc(input,mesh);
        else
            break;
        end
    end

elseif input.time_lapse_flag==1
 for ip_cnt=1:input.ip_num   
    for itr=1:input.itn
       [fem,mesh]=d4_prepare_data(itr,input,mesh,fem);
       [ex,fem,mesh]=rms_4d(itr,ip_cnt,input,mesh,fem);
       final=info_save(0,itr,ip_cnt,input,mesh,fem,final);
       if (ex==1 || ex==2) ;final.itr=itr-1;break; end
       [input]=update_lagran(itr,ip_cnt,1,input);
       [mesh]=kim_inversion2(input,mesh,fem);
       auto_contour_4d(itr,input,mesh,fem);
    end 
        if itr==input.itn+1; final.itr=itr;end
        if input.ip_flag==1
            [input,mesh]=ip_calc(input,mesh);
        else
            break;
        end   
 end
end

% Here save outputs...
info_save(3,itr,ip_cnt,input,mesh,fem,final);
save('inv_results.mat','final');

end