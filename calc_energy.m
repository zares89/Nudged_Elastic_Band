function [xq_int,yq_int,energy]=calc_energy(rep_mat,xq,yq,esurf_q)
% Linear interpolation to find energy on the path defined in rep_mat
for i=1:size(rep_mat,1)
    distx=[]; disty=[];
    for j=1:size(xq,2)
        distx(j)=abs(rep_mat(i,1)-xq(1,j));
    end
    [minabs,xq_id]=min(distx);
    for j=1:size(yq,1)
        disty(j)=abs(rep_mat(i,2)-yq(j,1));
    end
    [minabs,yq_id]=min(disty);
  %  xq_id=find(abs(rep_mat(i,1)-xq(1,:)) < 0.05);
   % yq_id=find(abs(rep_mat(i,2)-yq(:,1)) < 0.05);
    xq_int=xq_id; yq_int=yq_id;
    energy(i)=esurf_q(yq_int, xq_int);


end