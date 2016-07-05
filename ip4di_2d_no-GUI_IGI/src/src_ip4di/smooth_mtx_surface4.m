%
% TO DO: OUTPUT AND PLOT MATRICES CX, CY, CD1, CD2, MESH.M
%        (just to be sure of what they are)
%        THEN COMPUTE G'DG=(Cx Cy)D(Cx Cy)'
%        (first build D explicitly with a loop for getting u,v locally)
%

function mesh=smooth_mtx_surface4(input,mesh)

disp(' ')
disp('======================================')
disp('=     ENTER SMOOTH_MTX_SURFACE4      =')
disp('= define the model covariance matrix =')
disp('======================================')
disp(' ')

% init. covariance matrices
c=zeros(mesh.num_param,mesh.num_param);
cx=zeros(mesh.num_param,mesh.num_param);
cy=zeros(mesh.num_param,mesh.num_param);
cd1=zeros(mesh.num_param,mesh.num_param);
cd2=zeros(mesh.num_param,mesh.num_param);
mesh.ctc=zeros(mesh.num_param,mesh.num_param);
% tmp_x=union(mesh.tmp_param(:,1),mesh.tmp_param(:,1));
% tmp_y=union(mesh.tmp_param(:,2),mesh.tmp_param(:,2));

tmp_x=unique(mesh.param_x);
tmp_y=unique(mesh.param_y);


% build cx = 1st-order differential operator in the x-direction (x-gradient)
for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_x==current_x);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if ind~=length(tmp_x)
            if mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(ind+1) 
               cx(i,j)=1;  
               %cx(i,j)=sqrt(      (mesh.tmp_param(j,6)-mesh.tmp_param(j,5))/ ( mesh.tmp_param(j,1)-mesh.tmp_param(i,1)));
            end
        end
    end
end

for i=1:mesh.num_param
   cx(i,i)=-sum(cx(i,:));    
end

% covariance matrix
ctc1=cx'*cx;


% build cy
for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_y==current_y);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if ind~=length(tmp_y)
            if mesh.param_y(j)==tmp_y(ind+1) && mesh.param_x(j)==current_x 
               cy(i,j)=1;  
               %cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cy(i,i)=-sum(cy(i,:));    
end

% covariance matrix
ctc2=cy'*cy;           


% build cd1
for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    indy=find(tmp_y==current_y);
    indx=find(tmp_x==current_x);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if indy~=length(tmp_y) && indx~=length(tmp_x)
            if mesh.param_y(j)==tmp_y(indy+1) && mesh.param_x(j)==tmp_x(indx+1) 
               cd1(i,j)=1;  
               %cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cd1(i,i)=-sum(cd1(i,:));    
end

% covariance matrix
ctc3=cd1'*cd1;  %from top left to bottom right 


% build cd2
for i=1:mesh.num_param

    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    indx=find(tmp_x==current_x);
    indy=find(tmp_y==current_y);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if indy~=length(tmp_y) && indx~=1 
            if mesh.param_y(j)==tmp_y(indy+1) && mesh.param_x(j)==tmp_x(indx-1)
               cd2(i,j)=1;  
               %cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cd2(i,i)=-sum(cd2(i,:));    
end

% covariance matrix
ctc4=cd2'*cd2;   % from top right to bottom left 


% build total covariance matrix
if input.image_guidance==0
% use classic covariance matrix

  mesh.ctc=sparse(ctc1+ctc2+ctc3+ctc4);

  %%??
  %mesh.ctc=ctc1+ctc2;

elseif input.image_guidance==1
% use image-guided covariance via four-dir. smoothing (initcm1)

   [mesh]=initcm1(input,mesh,cx,cy,cd1,cd2);

elseif input.image_guidance==2
% use image-guidance covariance via local Laplacian (initcm2)

   [mesh]=initcm2(input,mesh,cx,cy);

elseif input.image_guidance==3
% read covariance matrix in input file

   [mesh.m2,mesh.m1]=size(mesh.map_param);

   fid=fopen(input.file_covariance,'r');
      mesh.ctc=fread(fid,[mesh.m1*mesh.m2,mesh.m1*mesh.m2],'real*4');
   fclose(fid);
end


%/* adjust weighting at the edges */

% if smooth_type==11
%     for i=1:num_param
%         s1=0;
%         for j=1:num_param
%             if(c(i,j)~=0) s1=s1+1; end
%         end
%         if(s1<=2) c(i,i)=1.15*c(i,j); end
%     end
% 
% end


%MK: This matrix has no meaning. I just calculate it so I can know the elements
%    that do not have zero. After this I can calcualte the S matrix (ACB, Kim)
%c=cx+cy;
c=mesh.ctc;

% /*Here I calculate the S matrix, which is one when C is non zero and zero otherwise*/
for i=1:mesh.num_param
    for j=1:mesh.num_param
        if (c(i,j)~=0)
           mesh.S(i,j)=1;
        else
           mesh.S(i,j)=0;
        end
    end
end


% Here create time related constrain matrix M
if input.time_lapse_flag==1

   %% JZ's matrix for third-order differential operator
   %%FL: dimensions are hard-coded, not sure how to use that...
   %mesh.M = spdiags([-ones(4125,1) 3*ones(4125,1) -3*ones(4125,1) ones(4125,1)],[0 375 750 1125],...
   %                  zeros(input.num_files*mesh.num_param,input.num_files*mesh.num_param));
   %mesh.M(3001:end,:)=0;

   %% Marios' matrix for first order differential operator
   mesh.M=eye(input.num_files*mesh.num_param,input.num_files*mesh.num_param);
   for i=1:input.num_files*mesh.num_param
       for j=1:input.num_files*mesh.num_param

           if j-i==mesh.num_param
              mesh.M(i,j)=-1;
           end

           if i==j && i>(input.num_files-1)*mesh.num_param && j>(input.num_files-1)*mesh.num_param
              mesh.M(i,j)=0;
           end
       end   %j
   end   %i

   %FL: here we keep derivative form to be able to apply the matrix to ACT in the future (not implemented yet),
   %    make sure to use the covariance form mesh.M'*mesh.M elsewhere (see rms_4d.m, kim_inversion2.m).

end   %end if time-lapse


% debug: plot matrices
if input.plot_covariance_matrices==1
   figure();
   imagesc(cx);
   title('Cx');
   colorbar();

   figure();
   imagesc(ctc1);
   title('C1=Cx''*Cx');
   colorbar();

   figure();
   imagesc(cy);
   title('Cy');
   colorbar();

   figure();
   imagesc(ctc2);
   title('C2=Cy''*Cy');
   colorbar();

   figure();
   imagesc(cd1);
   title('Cd1');
   colorbar();

   figure();
   imagesc(ctc3);
   title('C3=Cd1''*Cd1');
   colorbar();

   figure();
   imagesc(cd2);
   title('Cd2');
   colorbar();

   figure();
   imagesc(ctc4);
   title('C4=Cd2''*Cd2');
   colorbar();

   figure();
   imagesc(mesh.ctc);
   title('mesh.ctc');
   colorbar();

   if input.time_lapse_flag==1
      figure();
      imagesc(mesh.M);
      title('mesh.M');
      colorbar();

      figure();
      imagesc(mesh.M'*mesh.M);
      title('mesh.M'' x mesh.M');
      colorbar();
   end   %end if time-lapse

end

end   %end function mesh=smooth_mtx_surface4(input,mesh)

