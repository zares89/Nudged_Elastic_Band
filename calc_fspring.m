function [espring,fspring]=calc_fspring(rcoord,pcoord,rep_mat,k)

rep_mat_tmp=[rcoord;rep_mat;pcoord];
count=0;
for i=2:size(rep_mat_tmp,1)-1
    count=count+1;
    fspring(count,:)=k*(rep_mat_tmp(i+1,:)-rep_mat_tmp(i,:))-k*(rep_mat_tmp(i,:)-rep_mat_tmp(i-1,:));
    vec_forward= (rep_mat_tmp(i+1,:)-rep_mat_tmp(i,:));
    vec_backward= (rep_mat_tmp(i,:)-rep_mat_tmp(i-1,:));
    vec_tangent(count,:) = 0.5*(vec_forward+vec_backward);
    vec_tangent(count,:) = 1/norm(vec_tangent(count,:))*vec_tangent(count,:);
    fspring(count,:)=dot(fspring(count,:),vec_tangent(count,:))*vec_tangent(count,:);
    espring(count)=k/2*dot(rep_mat_tmp(i,:)-rep_mat_tmp(i-1,:),rep_mat_tmp(i,:)-rep_mat_tmp(i-1,:));
end
end