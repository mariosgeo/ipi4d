%
% FUNCTION EDIT_MODEL: define forward model from input file
%                      and eventually interpolate it on finite-element mesh
%
% Derived from function pushbutton1_Callback in M. Karaoulis' forward_modelling.m.
%
% Francois Lavoue', Colorado School of Mines, October 15, 2015

function [input,mesh]=edit_model(input,mesh)

disp(' ')
disp('---------------------------------------')
disp('         ENTER EDIT_MODEL              ')
disp(' read model parameters from input file ')
disp('---------------------------------------')
disp(' ')


faces2=[];
tmp1=1;

% first time call
if input.first_time_call==1

   % plot initial model (default values)
   if input.plot_model==1
      input.plot_options.label_title='DEFAULT MODEL';
      if input.dc_flag==1
         input.plot_options.cmplx_flag=0; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)
      elseif input.sip_flag==1
         input.plot_options.cmplx_flag=3; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)   %plot amplitude
         input.plot_options.cmplx_flag=4; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)   %plot phase
      end
   end

   if input.interp_model_flag==1

      disp('INTERPOLATE MODEL ON IRREGULAR GRID')
      disp(' ')

      % plot input model
      if input.plot_model==1
         input.plot_options.label_title='INPUT MODEL';
         input.plot_options.cmplx_flag=3; plot_model_on_input_grid(input,mesh,mesh.model_in);
         input.plot_options.cmplx_flag=4; plot_model_on_input_grid(input,mesh,mesh.model_in);
      end

      % reshape as rectangular matrix
      mesh.model_in=reshape(mesh.model_in,[input.n1,input.n2]);

      % input vectors
      vx1=[0:input.n2-1]*input.hstep;    %row
      vz1=[0:input.n1-1]'*input.hstep;   %column

      % output vectors
      vx2=mesh.XI(1,:);
      vz2=abs(mesh.YI(:,1));

      % complete initial vectors if output exceeds initial bounds
      if min(vx2)<min(vx1)
         vx1=[min(vx2) vx1];
         mesh.model_in=[mesh.model_in(:,1) mesh.model_in];
      end
      if max(vx2)>max(vx1)
         vx1=[vx1 max(vx2)];
         mesh.model_in=[mesh.model_in mesh.model_in(:,end)];
      end
      if min(vz2)<min(vz1)
         vz1=[min(vz2);vz1];
         mesh.model_in=[mesh.model_in(1,:);mesh.model_in];
      end
      if max(vz2)>max(vz1)
         vz1=[vz1;max(vz2)];
         mesh.model_in=[mesh.model_in;mesh.model_in(end,:)];
      end

      % interpolate (real and imag. parts separately for safety)
      tmp_real=interp2(vx1,vz1,real(mesh.model_in),vx2,vz2);
      tmp_imag=interp2(vx1,vz1,imag(mesh.model_in),vx2,vz2);
      mesh.res_param1=complex(tmp_real,tmp_imag);
      %size_res_param1=size(mesh.res_param1)
      %size_map_param =size(mesh.map_param)

      % recast model in single vector format
      % (be careful to non-conjugate transposition)
      mesh.res_param1=transpose(mesh.res_param1);
      mesh.res_param1=mesh.res_param1(:);
      %size_res_param1=size(mesh.res_param1)

      % plot
      if input.plot_model==1
         input.plot_options.label_title='INTERPOLATED MODEL';
         if input.dc_flag==1
            input.plot_options.cmplx_flag=0; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)
         elseif input.sip_flag==1
            input.plot_options.cmplx_flag=3; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)   %plot amplitude
            input.plot_options.cmplx_flag=4; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)   %plot phase
         end
      end

   end   %end interpolate model on irregular grid



   %
   % ROUND MODEL SUCH THAT IT CONTAINS ONLY TWO TYPES OF MEDIA (e.g. SAND 1, SAND 2)
   % (nearest interpolation in interp2 does not work...)
   if input.binary_model_flag==1

      disp('ROUND MODEL SUCH THAT IT CONTAINS ONLY TWO TYPES OF MEDIA (TYPE 1, TYPE 2)')

      val_ref1_r=max(real(mesh.res_param1));   %reference value for sand 1 (coarse, background)
      val_ref2_r=min(real(mesh.res_param1));   %    "       "    "   "   2 (fine, small-scale structures)

      val_ref1_i=max(imag(mesh.res_param1));
      val_ref2_i=min(imag(mesh.res_param1));

      for ipar=1:mesh.num_param

          % init. to very large values
          % (to detect bugs)
          rho_r=1e6;
          rho_i=1e6;

          % real part
          if abs(real(mesh.res_param1(ipar))-val_ref1_r)<=abs(real(mesh.res_param1(ipar))-val_ref2_r)
             rho_r=val_ref1_r;
          else
             rho_r=val_ref2_r;
          end

          % imaginary part
          if abs(imag(mesh.res_param1(ipar))-val_ref1_i)<=abs(imag(mesh.res_param1(ipar))-val_ref2_i)
             rho_i=val_ref1_i;
          else
             rho_i=val_ref2_i;
          end

          % update
          mesh.res_param1(ipar)=rho_r+1i*rho_i;

      end   %end loop over cells

   end   %end if binary model


   % re-define nb of parameters
   mesh.num_param=length(mesh.res_param1);

   % define background model
   input.bgr_res_param=mesh.res_param1;

   % plot
   if input.plot_model==1
      input.plot_options.label_title='FORWARD MODEL';
      if input.dc_flag==1
         input.plot_options.cmplx_flag=0; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)
      elseif input.sip_flag==1
         input.plot_options.cmplx_flag=3; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)   %plot amplitude
         input.plot_options.cmplx_flag=4; plot_model_on_forward_mesh(input,mesh,mesh.res_param1)   %plot phase
      end
   end

   input.first_time_call=0;

end   %end if first_time_call


%
% MANUAL EDIT OF THE MODEL (NOT USED)
%
input.edit_model=0;   %never edit model
if input.edit_model==1 then

[x y]=ginput(1);
% x=x-left(num_ext);
%bgr_res_param=[1:1:num_param];

if (x<min(mesh.tmp_param(:,3)) || x>max(mesh.tmp_param(:,4))) ||( y<-max(mesh.tmp_param(:,6)) ||  y>-min(mesh.tmp_param(:,5)))
   errordlg('Area selected if out of bounds','Selection Error');
else

   % Search which parameter user chose
   tmp_x=abs(mesh.param_x-x);

   [dummy,ind3]=min(tmp_x);
   x_need=mesh.param_x(ind3);
   [ind3]=find(mesh.param_x==x_need);

   tmp_y=abs( -y -mesh.param_y);

   [dummy,ind4]=min(tmp_y);
   y_need=mesh.param_y(ind4);
   [ind4]=find(mesh.param_y==y_need);

   % so we want the parameter withcoordinats (ind3,ind4);
   selected_param=intersect(ind3,ind4);
   % get new value for current parameter    
   if input.sip_flag==0
      current=inputdlg('Enter new resistivity value','Forward Model',1,{ num2str(input.bgr_res_param(selected_param))} );
      user_entry = str2double(current); 
   elseif input.sip_flag==1
      com_val{1}=num2str(abs(input.bgr_res_param(selected_param)));
      com_val{2}=num2str (  1000*atan2 (   imag(input.bgr_res_param(selected_param)),real(input.bgr_res_param(selected_param))    )   )  ;     
      current=inputdlg({'Enter new amplitude value','Enter new phase value'},'Forward Model',1, com_val );
      mag=str2double(current{1});
      phi=str2double(current{2});
      user_entry=complex (mag*(cos(phi/1000)) , mag*(sin(phi/1000)));
   end   %end sip_flag         
    
   %make some tests in case negative or nan
   if isnan(user_entry)
      errordlg('NaN','Not a Number','modal');
      uicontrol(hObject);
   elseif user_entry<0
      errordlg('Negative value','Negative Value','modal');
   elseif isempty(user_entry)==1
      % do nothing
   else
      disp('CHANGE PARAMETER')
      input.bgr_res_param(selected_param)=user_entry; 
   end

   tmp=mesh.faces{selected_param};
   tmp2=mesh.edge(tmp,:);
   tmp3=mesh.node(tmp2,:);
   tmp4=unique(tmp3,'rows');

   a=find(tmp4(:,2)==max(tmp4(:,2))); % up

   b=find(tmp4(:,2)==min(tmp4(:,2))); % up
   b=b(end:-1:1);

   c=find(tmp4(:,1)==min(tmp4(:,1))); % left
   cc=setdiff(c,a);
   c=setdiff(cc,b);

   d=find(tmp4(:,1)==max(tmp4(:,1))); % right
   dd=setdiff(d,a);
   d=setdiff(dd,b);

   faces=[a;d;b;c];

   patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',(abs( input.bgr_res_param(selected_param)))','edgecolor','k');
   plot(mesh.orig_probe_x,-mesh.orig_probe_z,'ko');
   %set(handles.pushbutton2,'Enable','On');    

end   %end if out of bounds

end   %end if edit model

