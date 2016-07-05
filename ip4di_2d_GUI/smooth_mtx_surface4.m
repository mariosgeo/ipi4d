function mesh=smooth_mtx_surface4(input,mesh)



c=zeros(mesh.num_param,mesh.num_param);
cx=zeros(mesh.num_param,mesh.num_param);
cy=zeros(mesh.num_param,mesh.num_param);

% tmp_x=union(mesh.tmp_param(:,1),mesh.tmp_param(:,1));
% tmp_y=union(mesh.tmp_param(:,2),mesh.tmp_param(:,2));


tmp_x=unique(mesh.param_x);
tmp_y=unique(mesh.param_y);



for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_x==current_x);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if ind~=length(tmp_x)
            if mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(ind+1) 
               cx(i,j)=1;  
%                cx(i,j)=sqrt(      (mesh.tmp_param(j,6)-mesh.tmp_param(j,5))/ ( mesh.tmp_param(j,1)-mesh.tmp_param(i,1)));
            end
        end
    end
end

for i=1:mesh.num_param
   cx(i,i)=-sum(cx(i,:));    
end


ctc1=cx'*cx;
        

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
            
            
 ctc2=cy'*cy;           

        
mesh.ctc=ctc1+ctc2;        
        

%ctc=eye(num_param,num_param);

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




% This matrix has no meaning. I just calculate so I can now the elemeents
% that do not have zero. After this I can calcualte the S matrix (ACB, Kim)
c=cx+cy;
c=mesh.ctc;
% /*Here i calculate the S matrix, which is one when C is non zero and zero otherwise*/
     for i=1:mesh.num_param
    
	  for j=1:mesh.num_param
           
            if (c(i,j)~=0)
                mesh.S(i,j)=1;
            else
                mesh.S(i,j)=0;
            end
      end
     end
    
     
     % Here create time related constarain matrix M
if input.time_lapse_flag==1
        mesh.M=eye(input.num_files*mesh.num_param,input.num_files*mesh.num_param);
        for i=1:input.num_files*mesh.num_param
            for j=1:input.num_files*mesh.num_param
        
            if j-i==mesh.num_param
                mesh.M(i,j)=-1;
            end

            if i==j && i>(input.num_files-1)*mesh.num_param && j>(input.num_files-1)*mesh.num_param
                mesh.M(i,j)=0;
            end
            end
        end

end


         
end        
        
