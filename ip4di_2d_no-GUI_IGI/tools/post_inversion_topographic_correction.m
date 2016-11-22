%
% This program correct an arbitrary input model for topography, using
% the 'mesh' structure stored during inversion. Useful for correcting
% re-interpolated models or training images for topography.
% 
% Francois Lavoue, 7 Oct. 2016

clear all
close all

%=====   USER PARAMETERS   =====%

% choose to use Matlab or Octave
% (does not work with Octave yet)
matlab_flag=1;   % 1->Matlab, else->Octave

% data format
data_format=3   % 2 -> (x,z)
                % 3 -> (x,z,v)

% choose inversion directory where to read topo information
suffix_out='_mesh-v2-15lay-807el_lagrn0.1_no-ACB_IGI-p01-plof4_force-st0'
dir=['results/inv3' suffix_out '/']

% file names
file_in=['results/inv3' suffix_out '/ellipses.dat']
file_out='ellipses_topo.dat'

file_in='data/training-image_profile1_Solfatara_n826x639_d0.7879x1.001896_xmin-3.018838_undersampled-x2_xzv.dat'
file_out='training-image_profile1_Solfatara_n826x639_d0.7879x1.001896_xmin-3.018838_undersampled-x2_topo_xzv.dat'

% structure file names
suffix_in='';
file_input_struct=[dir 'input_struct' suffix_in '.mat'];
file_mesh_struct =[dir 'mesh_struct'  suffix_in '.mat'];
file_final_struct =[dir 'final_struct'  suffix_in '.mat'];

% load structures from run dir.
if matlab_flag==1
   input=importdata(file_input_struct);
   mesh =importdata(file_mesh_struct);
   final=importdata(file_final_struct);
else   % Octave syntax
   load(file_input_struct);
   load(file_mesh_struct);
   load(file_final_struct);
end

%=====  END USER PARAMETERS   =====%


% read input model
data=load(file_in);
vMX=data(:,1);
vMZ=data(:,2);
if(data_format==3); model=data(:,3); end
ndata=length(vMX)

vx=sort(unique(vMX));
vz=sort(unique(vMZ));
nx=length(vx)
nz=length(vz)

xtopo=mesh.topo_nodes(:,1);
ztopo=mesh.topo_nodes(:,2);

% complete vector if necessary
if min(vx)<min(xtopo)
   xtopo = [ min(vx) ; xtopo];
   ztopo = [ ztopo(1); ztopo];
end
if max(vx)>max(xtopo)
   xtopo = [ xtopo ;  max(vx) ];
   ztopo = [ ztopo ;ztopo(end)];
end

% interpolate topo on input x-vector
ztopo_int=interp1(xtopo,ztopo,vx,'linear');

% correct for topography
for ix=1:nx
    ind=find(vMX==vx(ix));
    vMZ(ind)=vMZ(ind)-ztopo_int(ix);   %FL: WHY "-"???

    if data_format==3 && length(ind)~=nz
       error(sprintf('Found only %i points at x = %f (ix = %i), there should be nz = %i.',length(ind),vx(ix),ix,nz))
    end
end

% new vector of z-locations (does not match the size of model anymore)
vz=unique(vMZ);
nz2=length(vz)


% save model with topo correction
fid=fopen(file_out,'w');
for id=1:ndata

    if data_format==2
    % save (x,z)
       fprintf( fid,'%f %f \n',vMX(id),vMZ(id) );

    elseif data_format==3
    %save (x,z,v)
       fprintf( fid,'%f %f %f \n',vMX(id),vMZ(id),model(id) );

%        elseif input.cmplx_format==1
%        %save real and imaginary part of interpolated resistivity
%           fprintf( fid,'%f %f %f %f\n',...
%             vx(ix),vz(iz),real(model(iz,ix)),imag(model(iz,ix)) );
%        elseif input.cmplx_format==2
%        %save amplitude and phase of interpolated resistivity
%           fprintf( fid,'%f %f %f %f\n',...
%             vx(ix),vz(iz),abs(model(iz,ix)),1000*atan2(imag(model(iz,ix)),real(model(iz,ix))) );
    end   %iz
end   %ix
fclose(fid);

