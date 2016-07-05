%
% Function BUILD_DIFFERENTIAL_OPERATORS: build the 1st-order differential operators
%          needed by initcm1.m to define the 4-directional smoothing matrix
%          that serves as prior model covariance matrix.
%
% Author: Marios Karaoulis & Jieyi Zhou, Colorado School of Mines
% Version: July 2015
%
% NB: This code was initially in function smooth_mtx_surface4.m.
%     It has been moved here by F. Lavoue', September 2015,
%     to keep smooth_mtx_surface4.m simpler.

function [cx,cy,cd1,cd2]=build_differential_operators(mesh)

cx=zeros(mesh.num_param,mesh.num_param);   % 1st-order derivative in x-direction
cy=zeros(mesh.num_param,mesh.num_param);   %  "   "        "      in y-dir.
cd1=zeros(mesh.num_param,mesh.num_param);  %  "   "        "      in 1st diag. dir.
cd2=zeros(mesh.num_param,mesh.num_param);  %  "   "        "      in 2nd diag. dir.

% coordinates
tmp_x=unique(mesh.param_x);
tmp_y=unique(mesh.param_y);


% build cx
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
%                cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cy(i,i)=-sum(cy(i,:));
end


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
%                cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cd1(i,i)=-sum(cd1(i,:));
end


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
%                cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cd2(i,i)=-sum(cd2(i,:));
end
