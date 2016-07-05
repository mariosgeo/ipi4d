function mesh=smoothness_matrix(input,mesh)



% % cx=zeros(mesh.num_param,mesh.num_param);
% % cy=zeros(mesh.num_param,mesh.num_param);
% % cz=zeros(mesh.num_param,mesh.num_param);
% % 
% % x_centers=unique(mesh.param_i(:,1));
% % y_centers=unique(mesh.param_i(:,2));
% % z_centers=unique(mesh.param_i(:,3));
% % 
% % 
% % % weight factor
% % if mesh.scale_y>=1
% % 		 mesh.scaley=-1; mesh.scalex=-((1/mesh.scale_y)*(1/mesh.scale_y));
% % 		 mesh.scalez=-((mesh.scale_z/mesh.scale_y)*(mesh.scale_z/mesh.scale_y));
% % end
% % 
% % if mesh.scale_y<1
% % 		 mesh.scalex=-1; mesh.scaley=-((1/mesh.scale_x)*(1/mesh.scale_x));
% % % 		 mesh.scalez=-((1/mesh.scalex)*(1/mesh.scalex));
% %          mesh.scalez=-((mesh.scale_z/mesh.scale_x)*(mesh.scale_z/mesh.scale_x));
% % end
% %      
% % mesh.scalex=1;
% % mesh.scaley=1;
% % mesh.scalez=0.5;
% %       
% % % find cx
% % for i=1:mesh.num_param
% %     % current x_center is
% %     x_tmp=mesh.param_i(i,1);
% %     ind=find(x_centers==x_tmp);
% %     
% %     if ind~=length(x_centers) % exclude if we have the last center
% %         x_compare=x_centers(ind+1);
% %         diff=x_compare-x_tmp;
% %         cx(i,i)=mesh.scalex;
% %         for j=1:mesh.num_param
% %             % if same y and z center
% %             if (mesh.param_i(i,2)==mesh.param_i(j,2)) && (mesh.param_i(i,3)==mesh.param_i(j,3)) && mesh.param_i(j,1)-mesh.param_i(i,1)==diff
% %                 cx(i,j)=-mesh.scalex;
% %             end
% %         end
% %     end
% % end
% %     
% % ctc1=cx'*cx;
% % 
% % 
% % % find cy
% % for i=1:mesh.num_param
% %     % current y_center is
% %     y_tmp=mesh.param_i(i,2);
% %     ind=find(y_centers==y_tmp);
% %     
% %     if ind~=length(y_centers) % exclude if we have the last center
% %         y_compare=y_centers(ind+1);
% %         diff=y_compare-y_tmp;
% %         cy(i,i)=mesh.scaley;
% %         for j=1:mesh.num_param
% %             % if same y and z center
% %             if (mesh.param_i(i,1)==mesh.param_i(j,1)) && (mesh.param_i(i,3)==mesh.param_i(j,3)) && mesh.param_i(j,2)-mesh.param_i(i,2)==diff
% %                 cy(i,j)=-mesh.scaley;
% %             end
% %         end
% %     end
% % end
% % ctc2=cy'*cy;
% % 
% % 
% % 
% % 
% % % find cz
% % for i=1:mesh.num_param
% %     % current z_center is
% %     z_tmp=mesh.param_i(i,3);
% %     ind=find(z_centers==z_tmp);
% %     
% %     if ind~=length(z_centers) % exclude if we have the last center
% %         z_compare=z_centers(ind+1);
% %         diff=z_compare-z_tmp;
% %         cz(i,i)=mesh.scalez;
% %         for j=1:mesh.num_param
% %             % if same y and z center
% %             if (mesh.param_i(i,1)==mesh.param_i(j,1)) && (mesh.param_i(i,2)==mesh.param_i(j,2)) && mesh.param_i(j,3)-mesh.param_i(i,3)==diff
% %                 cz(i,j)=-mesh.scalez;
% %             end
% %         end
% %     end
% % end
% % ctc3=cz'*cz;
% % 
% % 
% % % c=cx+cy+cz;
% % 
% % 
% % mesh.ctc=ctc1+ctc2+ctc3;
% % mesh.c=cx+cy+cz;
% % % mesh.ctc=mesh.ctc'*mesh.ctc;
% % 
% % 
% % 
% % mesh.S=zeros(mesh.num_param,mesh.num_param);
% % % /*Here i calculate the S matrix, which is one when C is non zero and zero otherwise*/
% %      for i=1:mesh.num_param
% %     
% % 	  for j=1:mesh.num_param
% %            
% %             if (mesh.ctc(i,j)~=0)
% %                 mesh.S(i,j)=1;
% %             else
% %                 mesh.S(i,j)=0;
% %             end
% %       end
% %      end
% % mesh.cx=cx;
% % mesh.cy=cy;
% % mesh.cz=cz;
% % 
% % mesh.c=(mesh.cx+mesh.cy+mesh.cz);

%mesh.ctc=eye(mesh.num_param);


if input.time_lapse_flag==1
M=sparse([1:(input.num_files)*mesh.num_param],[1:(input.num_files)*mesh.num_param],[ones((input.num_files-1)*mesh.num_param,1) ;zeros(mesh.num_param,1)]);
M2=sparse([1:(input.num_files)*mesh.num_param],[1+mesh.num_param:(input.num_files)*mesh.num_param+mesh.num_param],[-1*ones((input.num_files-1)*mesh.num_param,1) ;zeros(mesh.num_param,1)]);
M2=M2(:,1:(input.num_files)*mesh.num_param);
mesh.M=M+M2;
end
% mesh.M=eye(input.num_files*mesh.num_param);
%         for i=1:input.num_files*mesh.num_param
%      % Here create time related constarain matrix M
% if input.time_lapse_flag==1
% %         mesh.M=sparse([1:(input.num_files)*mesh.num_param],[1:(input.num_files)*mesh.num_param],[ones((input.num_files-1)*mesh.num_param,1) zeros(mesh.num_param,1)]);
%         
%         for i=1:input.num_files*mesh.num_param
%             for j=1:input.num_files*mesh.num_param
%         
%             if j-i==mesh.num_param
%                 mesh.M(i,j)=-1;
%             end
% 
%             if i==j && i>(input.num_files-1)*mesh.num_param && j>(input.num_files-1)*mesh.num_param
%                 mesh.M(i,j)=0;
%             end
%             end
%         end
% 
% end
% 
%         end

clear cx cy cz ctc1 ctc2 ctc3
end