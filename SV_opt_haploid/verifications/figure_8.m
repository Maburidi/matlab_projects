clc
clear
close all 

% This script is to verify the 3D projections. It generates a couple of
% figures to show that 3D projection algorithm is making sense 


%%%%%%%%%  Projection of One 3D Non feasible Point %%%%%%%%%%%%%  

w=1;     
amp1=2;
phase2=0.5;
rand_3dpoints = amp1*rand(w,3)- phase2;
h=0;

for z=1:w

n=5000;     
nonfs_point = rand_3dpoints(z,:); 
[central_point] = Projections_3D(nonfs_point);

dis_3d_nonfs(z) = pdist2(central_point' ,nonfs_point); 

x_shift = central_point(1)*ones(n,1);
y_shift = central_point(2)*ones(n,1);
z_shift = central_point(3)*ones(n,1);

amp=0.3;
phase=0.15;
point = amp*rand(n,3)- phase;

point(:,1) = point(:,1) + x_shift; 
point(:,2) = point(:,2) + y_shift; 
point(:,3) = point(:,3) + z_shift; 

figure(1)
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3));hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;
str = ['+','o','*','s','d','^','p','h'];
%scatter3(point(:,1),point(:,2),point(:,3),'filled');hold on;
for ii=1:1
plot3(nonfs_point(1),nonfs_point(2),nonfs_point(3),str(ii),color="black",LineWidth=3); hold on 
plot3(central_point(1),central_point(2),central_point(3),str(ii),color="red",LineWidth=3);hold on
end 
xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
title("Random 3D Points")

xlim([-0.5 1.5])
ylim([-0.5 1.5])
zlim([-0.5 1.5])

hold off 

figure(2)
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3));hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;
str = ['+','o','*','s','d','^','p','h'];
scatter3(point(:,1),point(:,2),point(:,3),'filled');hold on;
for ii=1:1
plot3(nonfs_point(1),nonfs_point(2),nonfs_point(3),str(ii),color="black",LineWidth=3); hold on 
plot3(central_point(1),central_point(2),central_point(3),str(ii),color="red",LineWidth=3);hold on
end 
xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
title("Random 3D Points")
xlim([-0.5 1.5])
ylim([-0.5 1.5])
zlim([-0.5 1.5])

hold off 


% exclude non feasible points 
m=1; 
flag1 = 0; 
for i=1:n 
c = point(i,1); 
p = point(i,2);
gp = point(i,3); 

if 0 <= c && c <= 1 && 0 <= p && p <= 1 && 0<=gp && gp<=1 && c <= p && p<=gp   %Region inside the triangle
            point2(m,1) = c;
            point2(m,2) = p;
            point2(m,3) = gp;
            feas_po= point2(m,:);
            dis(m,z) = pdist2(nonfs_point,feas_po);
            m=m+1;
else
    continue
end
end

num_fspoi(z) = m-1; 
tst = dis(1:num_fspoi(z),z) < dis_3d_nonfs(z) ;

if sum(tst) > 0
flag1 = flag1+1; 
disp("===== Failed ======")
disp(z)
end 

end 


figure(3)
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3),'.');hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;
str = ['+','o','*','s','d','^','p','h'];
%scatter3(all_poi(:,1),all_poi(:,2),all_poi(:,3),"filled",[0 0 0 ; 1 1 1; 4 4 4]');
scatter3(point2(:,1),point2(:,2),point2(:,3),'filled');hold on;

for ii=1:1
plot3(nonfs_point(1),nonfs_point(2),nonfs_point(3),str(ii),color="black",LineWidth=2); hold on 
plot3(central_point(1),central_point(2),central_point(3),str(ii),color="red",LineWidth=2);hold on
end 

xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
xlim([-0.5 1.5])
ylim([-0.5 1.5])
zlim([-0.5 1.5])

hold off 

%%% Plot results: 
figure(4)

for u=1:w
x_ = u*ones(num_fspoi(u),1) ;
y_= dis(1:num_fspoi(u),u) ; 
scatter(x_, y_ ,  'filled','blue');  hold on                       
end 
scatter(1:w, dis_3d_nonfs, 'filled','red');hold off                         
xlabel('$Point$','FontSize',30,Interpreter="latex");
ylabel('Distance from Origional Point','FontSize',20);






%%


%%%%%%%%%  Projection of 100 3D Non-feasible Point %%%%%%%%%%%%%  
% chnage w to increase the number of the non feasible points 

clc
clear
close all 
%%%%% Genertae Random 3D Points 

w=100;     
amp1=2;
phase2=0.5;
rand_3dpoints = amp1*rand(w,3)- phase2;
h=0;

for z=1:w

n=5000;     
nonfs_point = rand_3dpoints(z,:); 
[central_point] = Projections_3D(nonfs_point);

dis_3d_nonfs(z) = pdist2(central_point' ,nonfs_point); 

x_shift = central_point(1)*ones(n,1);
y_shift = central_point(2)*ones(n,1);
z_shift = central_point(3)*ones(n,1);

amp=0.3;
phase=0.15;
point = amp*rand(n,3)- phase;

point(:,1) = point(:,1) + x_shift; 
point(:,2) = point(:,2) + y_shift; 
point(:,3) = point(:,3) + z_shift; 


% exclude non feasible points 
m=1; 
flag1 = 0; 
for i=1:n 
c = point(i,1); 
p = point(i,2);
gp = point(i,3); 

if 0 <= c && c <= 1 && 0 <= p && p <= 1 && 0<=gp && gp<=1 && c <= p && p<=gp   %Region inside the triangle
            point2(m,1) = c;
            point2(m,2) = p;
            point2(m,3) = gp;
            feas_po= point2(m,:);
            dis(m,z) = pdist2(nonfs_point,feas_po);
            m=m+1;
else
    continue
end
end

num_fspoi(z) = m-1; 
tst = dis(1:num_fspoi(z),z) < dis_3d_nonfs(z) ;

if sum(tst) > 0
flag1 = flag1+1; 
disp("===== Failed ======")
disp(z)
end 

end 

%%% Plot results: 
figure(1)

for u=1:w
x_ = u*ones(num_fspoi(u),1) ;
y_= dis(1:num_fspoi(u),u) ; 
scatter(x_, y_ ,  'filled','blue');  hold on                       
end 
scatter(1:w, dis_3d_nonfs, 'filled','red');hold off                         
xlabel('Point Number','FontSize',20);
ylabel('Distance from Non-Feasible Point','FontSize',20);


%%
function [f_subs]= Projections_3D(f_vec)

N = length(f_vec);
n = N/3; 

c_vec = f_vec(1:n);
p_vec = f_vec(n+1:2*n);
gp_vec = f_vec(2*n+1:3*n);

%Initialize vectors r, s, and t
r = zeros(n,1);
s = zeros(n,1);
t = zeros(n,1);
%region=zeros(m,1);


% This ....
for i=1:n
    c=c_vec(i);
    p=p_vec(i);
    gp=gp_vec(i);


%--------------------------------------------------
% Inside Region
%--------------------------------------------------
    if   0 <= c && c <= 1 && 0 <= p && p <= 1 && 0<=gp && gp<=1 && c <= p && p<=gp 
            r(i)=c;
            s(i)=p;
            t(i)=gp;


    elseif   -gp + p <= 0 && -c <= 0 && c - p <= 0 && -1 + gp <= 0
            r(i)=c;
            s(i)=p;
            t(i)=gp;
    elseif  sqrt(2)*c >= 0 && sqrt(2) - gp/sqrt(2) - p/sqrt(2) >= 0 && ...
            -sqrt(2) - sqrt(2)*(-1 + c) + gp/sqrt(2) + p/sqrt(2) >= 0 && -gp + p > 0
            r(i)=c;
            s(i)=c + (1/2)*(-2*c + gp + p);
            t(i)=c + (1/2)*(-2*c + gp + p); 

     elseif sqrt(2)*(gp/sqrt(2) - p/sqrt(2)) >= 0 && ...
            -((-(gp/sqrt(2)) - p/sqrt(2))/sqrt(2)) - (gp/sqrt(2) - p/sqrt(2))/sqrt(2) >= 0 && ...
            (1/sqrt(2) - gp/sqrt(2) - p/sqrt(2))/sqrt(2) + (1/sqrt(2) - sqrt(2))*(-(1/sqrt(2)) + gp/sqrt(2) - p/sqrt(2)) >= 0 && ...
            -c > 0
            r(i)=0;
            s(i)=p;
            t(i)=gp;  

     elseif c/sqrt(2) + p/sqrt(2) >= 0 && -(c/sqrt(2)) + sqrt(2)*gp - p/sqrt(2) >= 0 && ...
            sqrt(2)*(1 - gp) >= 0 && c - p > 0 
            r(i)=(c+p)/2;                       
            s(i)=(c+p)/2;                       
            t(i)= gp;                 

     elseif 1 - p >= 0 && c >= 0 && -c + p >= 0 && -1 + gp > 0 
            r(i)=c;
            s(i)=p;
            t(i)=1;
    elseif 2*c < 0 && -gp + p > 0 && gp + p >= 0 && 2 - gp - p >= 0
            r(i)=0;
            s(i)=gp/2 + p/2;
            t(i)=gp/2 + p/2; 


    elseif p < 0 && -c - p > 0 && gp >= 0 && 1 - gp >= 0
            r(i)=0;
            s(i)=0;
            t(i)= gp;

    elseif 2*c - gp - p > 0 && -c + 2*gp - p < 0 && c + gp + p >= 0 && ... 
           3 - c - gp - p >= 0
            r(i)=c/3 + gp/3 + p/3;
            s(i)=c/3 + gp/3 + p/3;
            t(i)= c/3 + gp/3 + p/3; 

    elseif  -1 + gp > 0 && c < 0 && 1 - p >= 0 && p >= 0  
            r(i)=0;
            s(i)=p;
            t(i)=1; 

    elseif 2*(-1 + gp) > 0 && -c + p < 0 && c + p >= 0 && 2 - c - p >= 0
            r(i)=c/2 + p/2;
            s(i)=c/2 + p/2;
            t(i)=1; 

    elseif 2 - gp - p < 0 && -1 + p > 0 && c >= 0 && 1 - c >= 0
            r(i)=c;
            s(i)=1;
            t(i)=1; 

    elseif gp + p < 0 && gp < 0 && c + gp + p < 0
           r(i)=0;
           s(i)=0;
           t(i)=0; 
            
    elseif  2 - gp - p < 0 && 1 - p < 0 && c < 0
            r(i)=0;
            s(i)=1;
            t(i)=1;  
    elseif  1 - gp < 0 && p < 0 && c + p < 0
            r(i)=0;
            s(i)=0;
            t(i)=1;  
    else    
            r(i)= 1;    
            s(i)= 1;   
            t(i)= 1;   
    end 
end 
f_subs=[r;s;t];
end 






