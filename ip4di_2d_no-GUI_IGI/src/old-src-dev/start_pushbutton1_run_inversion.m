% --- Executes on button press in pushbutton1.
%function pushbutton1_Callback(hObject, eventdata, handles)
function start_pushbutton1_run_inversion
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mesh input tim

if input.time_lapse_flag==1
    figure;
    if input.dc_flag==1
        for i=1:input.num_files
            tim(i)=subplot(1,input.num_files,i);
        end
    elseif input.ip_flag==1 || input.sip_flag==1
        for i=1:2*input.num_files
            tim(i)=subplot(2,input.num_files,i);
        end
    end
end

% INVERSION
main(input,mesh)

