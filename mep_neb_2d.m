function mep_neb_2d()
%% Deriving minimum energy path from nudged elastic band algorithm (Henkelman and Jonsson, J Chem Phys, 113, 9978-9985 (2000)
%% Read free energy data
txtfile=strcat(['freefile-all-matlab-revised-3.dat']);
mydata1=importdata(txtfile);
x_dist=mydata1(:,1); % Oh distance from the surface
y_dist=mydata1(:,2); % Surface magnesium coordination number
freeval=mydata1(:,3); % Free energy (kcal/mol)

%% Store free energy values in the xy plane in the matrix 'freeval'
free_mat=[];
count=0;
free_cell={};
x_dist_old=x_dist(1);
countx=0; % x index
ynum=[];
x_list=[];
y_list=[];
while count < length(x_dist)
    
    count=count+1;
   
    if x_dist(count) ~= x_dist_old
        countx=countx+1;
        free_cell{countx}=free_mat;
        x_list=[x_list;x_dist_old];
        x_dist_old=x_dist(count);
        y_num(countx)=length(free_mat);
        free_mat=[];
        
    elseif count==length(x_dist)
        countx=countx+1;
        free_cell{countx}=[free_mat;freeval(count)];
         y_num(countx)=length(free_mat)+1;
         x_list=[x_list;x_dist_old];
    end
    free_mat=[free_mat;freeval(count)];
    
end
y_list=y_dist(1:y_num(1));
freeval=zeros(countx,y_num(1));
for i=1:countx
    for j=1:y_num(1)
        if free_cell{i}(j) > 25 || j >= 45 
          freeval(i,j)=nan;
        elseif i < 15 && j > 31
            freeval(i,j)=nan;
        else
        freeval(i,j)=free_cell{i}(j)-2.101;
        end
    end
end
%% Plot the free energy landscape
figure (15)
 newpoints_x =100;
 newpoints_y = 100;
  cid_start=4;
  xcid_start=2;
  y_listp=y_list(cid_start:end);
  x_listp=x_list(xcid_start:end);
[xq,yq] = meshgrid(...
            linspace(min(min(x_listp,[],2)),max(max(x_listp,[],2)),newpoints_x ),...
            linspace(min(min(y_listp,[],1)),max(max(y_listp,[],1)),newpoints_y )...
          );
esurf_q = interp2(x_listp,y_listp,freeval(xcid_start:end,cid_start:end)',xq,yq,'linear');
 cMap=jet(256); %set the colomap using the "jet" scale
[c,h]=contourf(xq-7.55-0.096,yq,esurf_q)
set(h, 'edgecolor','none');
colormap(cMap);
hold on
xlabel('Z (\AA)')
ylabel('Mg^{2+} Coordination Number')
title('Free energy kcal/mol')
colorbar
 saveas(figure(15),'2_20_free.fig');
 %figure (16)
 %surf(x_list,y_list,freeval');
 %saveas(figure(16),'2_20_free3D_2CV.fig');
 
 
 %% Calculate Minimum Energy Path (MEP) from a simple nudged elastic band algorithm (Henkelman and Jonsson, J Chem Phys, 113, 9978-9985 (2000).
 % make a mesh in x and y direction and interpolate 

%figure (17)
%surf(xq,yq,esurf_q);
% Reactant and Product coordinates
rcoord=[7.65 2.35];pcoord=[10.25 4.05];

% Linear interpolation between initial and final states with nreplica (number
% of replicas)
nreplica=21;
x_grid=linspace(rcoord(1),pcoord(1),nreplica);
x_grid=x_grid(2:length(x_grid)-1);
y_grid=linspace(rcoord(2),pcoord(2),nreplica);
y_grid=y_grid(2:length(y_grid)-1);
rep_mat=[x_grid' y_grid'];
%plot(rep_mat(:,1),rep_mat(:,2));
%hold on;
% NEB and Minimization parameters
k=2000;imax=10000;sigma0=0.300;jmax=1000;eps_s=0.001;


[xq_int,yq_int,intenergy]=calc_energy(rep_mat,xq,yq,esurf_q);% Find energy values on the linear path defined above in rep_mat
[espring,fspring]=calc_fspring(rcoord,pcoord,rep_mat,k); % Find spring force and energy (fspring, fenergy)
fland=calc_fland(rep_mat,xq,yq,esurf_q,rcoord,pcoord); % Find true force (fland)
fnudged=fland+fspring; % Total force = fspring + fland
pot_old=0;
pot_new=0;
tol=10;
direc=fnudged;
for i=1:size(fnudged,1)
    pot_new=pot_new+intenergy(i)+espring(i);
end
rep_mat_new=rep_mat;
% In the itartions below, replicas are moved according to fnudged force until optimization is achieved 
while i <= imax && tol >= 0.0001
    jj=0;
    pot_old=pot_new;
    rep_mat_old=rep_mat_new;
    gradnorm=0;
    %gradmat= calc_grad(fnudged,rep_mat_new);
    for j=1:size(rep_mat_new,1)
        gradnorm=gradnorm+dot(fnudged(j,:),fnudged(j,:));
    end
    gradnorm=sqrt(gradnorm);
    
    %
    alpha=-sigma0;
    fold=fnudged;
    
    for j=1:size(rep_mat_new,1)
        rep_mat_new(j,:)=rep_mat_old(j,:)+sigma0*direc(j,:)/gradnorm;
    end
    [espring,fspring]=calc_fspring(rcoord,pcoord,rep_mat_new,k);
    fland=calc_fland(rep_mat_new,xq,yq,esurf_q,rcoord,pcoord);
    fnudged=fland+fspring;
    %
    eta_prev=0;
    for j=1:size(rep_mat_new,1)
        eta_prev=eta_prev+dot(fnudged(j,:),direc(j,:));
    end
    %
    while jj <= jmax && alpha*alpha*gradnorm*gradnorm >= eps_s*eps_s
        eta=0;
        for j=1:size(rep_mat_new,1)
            eta=eta+dot(fold(j,:),direc(j,:));
        end
        alpha=alpha*eta/(eta_prev-eta);
        for j=1:size(rep_mat_new,1)
            rep_mat_new(j,:)=rep_mat_old(j,:) + alpha*direc(j,:)/gradnorm;
            rep_mat_old(j,:)= rep_mat_new(j,:);
        end
        [espring,fspring]=calc_fspring(rcoord,pcoord,rep_mat_new,k);
        fland=calc_fland(rep_mat_new,xq,yq,esurf_q,rcoord,pcoord);
        fnudged=fland+fspring;
        [xq_int,yq_int,intenergy]=calc_energy(rep_mat_new,xq,yq,esurf_q);
        eta_prev=eta;
        jj=jj+1;
        pot_new=0;
        fold=fnudged;
        for j=1:size(rep_mat_new,1)
           % fold(j,:)=fnudged(j,:);
            pot_new=pot_new+intenergy(j)+espring(j);
        end
        
        %for j=1:size(rep_mat_new,1)
       %     direc(j,:)=fnudged(j,:);
      %  end
        
    end
    direc=fnudged;
    tol=abs(pot_new-pot_old)
    i=i+1;
    %plot(rep_mat_new(:,1),rep_mat_new(:,2));
%hold on;
end
rep_mat_path1=rep_mat_new;
rcoord=[7.65 2.35];pcoord=[10.25 4.05];
x_grid=linspace(rcoord(1),pcoord(1),33);
x_grid=x_grid(2:length(x_grid)-1);
y_grid=linspace(rcoord(2),pcoord(2),33);
y_grid=y_grid(2:length(y_grid)-1);
rep_mat=[x_grid' y_grid'];
%plot(rep_mat(:,1),rep_mat(:,2));
%hold on;

rep_mat_all=[rep_mat_path1];
rep_mat_allt=rep_mat_all;
rep_mat_allt(:,1)=rep_mat_allt(:,1)-7.55-0.096
rep_mat_allt=[7.65-7.55-0.096 2.35;rep_mat_allt;10.25-7.55-0.096 4.05];
save rep.mat rep_mat_all
 plot(rep_mat_allt(:,1),rep_mat_allt(:,2),'-ks',...
    'LineWidth',2,...
    'Marker','o','MarkerSize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');


rcoord=[7.65 2.35];pcoord=[10.25 4.05];
 [fe,xval,yval]=calc_energy_all(rep_mat_all,xq,yq,esurf_q,rcoord,pcoord);
 num_replicas=21;
 figure (123)
 plot(1:num_replicas,fe);
 %hold on
 %plot(1:num_replicas,fe0);
 xlabel('\Replica}')
ylabel('\Free Energy (kcal/mol)')
title('MEP')
colorbar
 saveas(figure(123),'2_20_MEP1-k400-r43.fig');
 file=fopen('x_energy.txt','w');
for i=1:length(xval)
    fprintf(file,'%4.8f %4.8f\n',i,fe(i));
end
fclose(file);
 file=fopen('x_coord.txt','w');
for i=1:length(xval)
    fprintf(file,'%4.8f %4.8f\n',i,yval(i));
end
fclose(file);

%figure (233)


%f_nudged=  
%  xstart=[9.85 2.55];x_end=[2.05 1.55];
% xstart_id=find(x_list==xstart(1));
%  ystart_id=find(y_list==xstart(2));
%  xend_id=find(x_list==x_end(1));
%  yend_id=find(y_list==x_end(2));
%  num_x=abs(xend_id-xstart_id);%num_x=70;
%  num_y=abs(yend_id-ystart_id);
%  dec_max=0;
%  allpath=zeros(factorial(num_x+num_y)/(factorial(num_x)*factorial(num_y)),num_x+num_y);
%  for i=1:num_x
%      for k=0:1
%          
%          for j=1:num_y
%              for l=0:1
%                  allpath((i-1)*2*num_x+k+(j-1)*num_y*2+l+1,i+j-1)=[xstart_id xstart_id+k
%              end
%          end
%      end
%  end
%  for i=1:num_y
%      dec_max=sym(dec_max+2^(i-1+num_x));
%  end
%  dec_min=0;
%   for i=1:num_y
%      dec_min=dec_min+2^(i-1);
%   end
%   arr_bin=de2bi(double(dec_max));
%   mat_bin=zeros(1,length(arr_bin));
%   count=0;
%   i=dec_max;
%   while i >= dec_min
% % for i=dec_max:-1:dec_min
%      count=count+1;
%      arr_bin=de2bi(double(i));
%      width=length(arr_bin);
%      mat_bin(count,length(arr_bin)-width+1:length(arr_bin))=arr_bin;
%      i=i-1;
%  end
 
end