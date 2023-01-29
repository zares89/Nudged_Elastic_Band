function fland=calc_fland(rep_mat,xq,yq,esurf_q,rcoord,pcoord)
% Calculation of the true force on the band between rcoord and pcoord defined
% in rep_mat
rep_mat_tmp=[rcoord;rep_mat;pcoord];
count=0;

for i=2:size(rep_mat_tmp,1)-1
    count=count+1;
    distx=[]; disty=[];
    for j=1:size(xq,2)
        distx(j)=abs(rep_mat_tmp(i,1)-xq(1,j));
    end
    [minabs,xq_id]=min(distx);
    for j=1:size(yq,1)
        disty(j)=abs(rep_mat_tmp(i,2)-yq(j,1));
    end
    [minabs,yq_id]=min(disty);
    xq_int=xq_id(1); yq_int=yq_id(1);
    energy(count)=esurf_q(yq_int, xq_int);
    % Calculate gradient based on forward-backward finite difference
    xq_forward=xq_int+1;xq_backward=xq_int-1;
    gradxf=(esurf_q(yq_int,xq_forward)-energy(count))/(xq_forward-xq_int);
    gradxb=(-esurf_q(yq_int,xq_backward)+energy(count))/(-xq_backward+xq_int);
    gradx=0.5*(gradxf+gradxb);
    yq_forward=yq_int+1;yq_backward=yq_int-1;
    gradyf=(esurf_q(yq_forward,xq_int)-energy(count))/(yq_forward-yq_int);
    gradyb=(-esurf_q(yq_backward,xq_int)+energy(count))/(-yq_backward+yq_int);
    grady=0.5*(gradyf+gradyb);
    gradvec(count,:)=[gradx grady];
    % Tangential vector along the path
    tangf=rep_mat_tmp(i+1,:)-rep_mat_tmp(i,:);
    tangb=rep_mat_tmp(i,:)-rep_mat_tmp(i-1,:);
    tang(count,:)=0.5*(tangf+tangb);
    tang(count,:)=1/norm(tang(count,:))*tang(count,:);
    % F(Ri) force perpendicular to the elastic band
    fland(count,:)=-gradvec(count,:)+dot(gradvec(count,:),tang(count,:))*tang(count,:);
end

end