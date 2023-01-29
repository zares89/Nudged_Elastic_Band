function [fe,xval,yval]=calc_energy_all(rep_mat,xq,yq,esurf_q,rcoord,pcoord)

rep_mat_tmp=[rcoord;rep_mat;pcoord];

for i=1:size(rep_mat_tmp,1)
    distx=[]; disty=[];
    for j=1:size(xq,2)
        distx(j)=abs(rep_mat_tmp(i,1)-xq(1,j));
    end
    [minabs,xq_id]=min(distx);
    for j=1:size(yq,1)
        disty(j)=abs(rep_mat_tmp(i,2)-yq(j,1));
    end
    [minabs,yq_id]=min(disty);
  %  xq_id=find(abs(rep_mat(i,1)-xq(1,:)) < 0.05);
   % yq_id=find(abs(rep_mat(i,2)-yq(:,1)) < 0.05);
    xq_int(i)=xq_id; yq_int(i)=yq_id;
    fe(i)=esurf_q(yq_int(i), xq_int(i));
    xval(i)=rep_mat_tmp( i,1);
    yval(i)=rep_mat_tmp( i,2);
end