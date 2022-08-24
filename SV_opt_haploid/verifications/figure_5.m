% 2D Projections verification 
%%%% this script is to implement verifications and comparisons of 2D vs 3D
%%%% projection approaches. 
%%%% It generate 7 figures 

clc
clear
close all 
%%%%% Genertae Random 3D Points 
n=10000;     
amp=2;
phase=0.5;
point = amp*rand(n,3)- phase;

figure(1)

corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3),'.');hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;

scatter3(point(:,1),point(:,2),point(:,3),'filled');hold off;
xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
title("Random 3D Points")

xlim([-0.5 1.5])
ylim([-0.5 1.5])
zlim([-0.5 1.5])

hold off 

%% Project the points - on 2D and 3D approaches/compute the distances   
d=0;
m=0;
for i=1:n
    [f_3dproj] = Projections_3D(point(i,:));
    [f_2dproj] = Alt2D_Proj_II(point(i,:),30);
    points_2d(i,:)= f_2dproj;
    points_3d(i,:)= f_3dproj;
    dist_from_2d(i) = pdist2(f_2dproj' , point(i,:));
    dist_from_3d(i) = pdist2(f_3dproj' , point(i,:));
    dis_2d_3d(i) = pdist2(f_2dproj',f_3dproj'); 
    dis2d_dis3d(i) = dist_from_2d(i) - dist_from_3d(i); 

    if abs(dis2d_dis3d(i)) > 0.00001     
        m=m+1;                           
        unsu_pots_3d(m,:) = f_3dproj;    
        unsu_pots_2d(m,:) = f_2dproj;    
        unus_dist(m) = dis2d_dis3d(i);    
    end 

    % if (dis2d_dis3d(i) < -10^(-15))
     %   d=d+1;
      %  unsu_pots_3d(d,:) = f_3dproj;
       % unsu_pots_2d(d,:) = f_2dproj;
        % unus_dist(d) = dis2d_dis3d(i); 
         % end 
end 
%unus_dist =unus_dist';





%%  Plot Projections 
figure (2)                 
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3));hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;
%for zz =1:ii
%scatter3(pots_3d_dif(:,1),pots_3d_dif(:,2),pots_3d_dif(:,3),'filled','MarkerFaceColor',[0.75 .75 .75]);hold on;
%scatter3(pots_2d_dif(:,1),pots_2d_dif(:,2),pots_2d_dif(:,3),'filled','MarkerFaceColor',[0 .1 .1]);hold on;
scatter3(points_3d(:,1),points_3d(:,2),points_3d(:,3),'filled','MarkerFaceColor',[0 .1 .1]);hold on;

xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
title("3 D Projections")
hold off 

figure(3)
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3));hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;
scatter3(points_2d(:,1),points_2d(:,2),points_2d(:,3),'filled','MarkerFaceColor',[0 .1 .1]);hold on;

%end 
xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
title("2 D Projections")
hold off


%% Compare some points 
figure(4)
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3),'.');hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=1, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=1, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=1, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=1, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=1, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=1, Color="red"); hold on;
str = ['+','o','*','s','d','^','p','h'];
for ii=1:1%length(str) 
%disp(point(ii,:))
%disp(points_2d(ii,:))
%disp(points_3d(ii,:))
%disp(dist_from_3d(ii))
%disp(dist_from_2d(ii))
%disp("-----------------------")
plot3(point(ii,1),point(ii,2),point(ii,3),str(ii),color="black",LineWidth=2); hold on 
plot3(points_2d(ii,1),points_2d(ii,2),points_2d(ii,3),str(ii),color="blue",LineWidth=2);hold on
plot3(points_3d(ii,1),points_3d(ii,2),points_3d(ii,3),str(ii),color="red",LineWidth=2);hold on
end 

%scatter3(all_poi(:,1),all_poi(:,2),all_poi(:,3),"filled",[0 0 0 ; 1 1 1; 4 4 4]');
xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
xlim([-0.5 1.5])
ylim([-0.5 1.5])
zlim([-0.5 1.5])

%legend("point","2D","3D")
hold off

%%  Plot Distance 1
% between the 2D and origional point, and 3D and origional point projections  

figure(5)
num_poi=300;
stem( dist_from_2d(1:num_poi), 'filled','red');    hold on;
stem( dist_from_3d(1:num_poi), 'filled','blue')    

xlabel('Point Number','FontSize',30,Interpreter="latex");
ylabel('Euclidean Distance to Original Point','FontSize',20, Interpreter="latex");
legend("2D","3D") 
ylim([-0.01 2.5]) 
hold off 

%% Plot Distance 2
% between point projected by 2D and point projected by 3D   
figure(6)
scatter(1:n, dis_2d_3d, 'filled','red')                         
xlabel('Point Number','FontSize',30,Interpreter="latex");
ylabel('Euclidean Distance','FontSize',30, Interpreter="latex");
title("Distance between point projected by 2D and same point projected by 3D")
dim1 = [.2 .5 .3 .37];
str1 = strcat('Mathings Percent= ',num2str(length(find( dis_2d_3d < 10^(-10) ))/n));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');

ylim([-0.3 2])
hold off 
%figure(7) 
%histfit(dis_2d_3d)

%% Plot dist2d - dist3d
figure(7)
scatter(1:n, dis2d_dis3d, 'filled','red')
xlabel('Point Number','FontSize',30,Interpreter="latex");
ylabel('Difference between distances','FontSize',30, Interpreter="latex");
%title("Distance between point projected by 2D and same point projected by 3D")
%dim1 = [.2 .5 .3 .37];
%str1 = strcat('Mathings Percent= ',num2str(length(find( dis_2d_3d ==0))/n));
%annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
ylim([-0.3 2])
hold off 

%%

figure (8)                 
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3));hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=3, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=3, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=3, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=3, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=3, Color="red"); hold on;
%for zz =1:ii
%scatter3(pots_3d_dif(:,1),pots_3d_dif(:,2),pots_3d_dif(:,3),'filled','MarkerFaceColor',[0.75 .75 .75]);hold on;
%scatter3(pots_2d_dif(:,1),pots_2d_dif(:,2),pots_2d_dif(:,3),'filled','MarkerFaceColor',[0 .1 .1]);hold on;
scatter3(unsu_pots_3d(:,1),unsu_pots_3d(:,2),unsu_pots_3d(:,3),'filled','MarkerFaceColor',[0 .1 .1]);hold on;

xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
title("3 D Projections")
hold off 








%%

function [a_new,b_new] = Projec_2D(a,b,c) 

   if  0<=a && a <= b && 0 <= b  && b <= c            % Check if the point is interior
        a_new = a; 
        b_new = b; 
   elseif a <= 0 && 0 <= b  && b <= c   %Region (1)
        a_new = 0;
        b_new = b;
   elseif 0 <= a && a <= c && b >= c    %Region (2)
        a_new = a;
        b_new = c;
   elseif a <= 0 && b >= c               %Region (3)
        a_new = 0;
        b_new = c;
    elseif a <= -b && b <= 0            %Region (4)
        a_new = 0;
        b_new = 0;
    elseif a >= c && b >= (2*c - a)       %Region (5)
        a_new = c;
        b_new = c;
    else                               %Regions (6)
        a_new = (a + b)/2;
        b_new = (a + b)/2;
    end
end 

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

function [ f_subs] = Alt2D_Proj_II(f_vec,maxit) 
thresh=0.000000000001; 
N = length(f_vec);
n = N/3 ;
c_vec = f_vec(1:n);
p_vec = f_vec(n+1:2*n);
gp_vec = f_vec(2*n+1:3*n);

r = zeros(n,1);
s = zeros(n,1);
t = zeros(n,1);
for i = 1:n
    c = c_vec(i);
    p = p_vec(i);
    gp = gp_vec(i);

    if 0 <= c && c <= 1 && 0 <= p && p <= 1 && 0<=gp && gp<=1 && c <= p && p<=gp   %Region inside the triangle
               r(i) = c;
               s(i) = p;
               t(i) = gp;
              continue;
    else     
        j =1; 
        while (j < maxit)
             % if j==1
             % zppos = sort([0,gp,1]);
             % gp = zppos(2);
             % end
             %----------------------------------------------------------------------
               [c_new, p_new] = Projec_2D(c, p, 1);
               [c_new, gp_new] = Projec_2D(c_new, gp, 1); 
               [p_new, gp_new] = Projec_2D(p_new, gp_new, 1);
             % disp(j)
                j = j+1;
             %----------------------------------------------------------------------
              if abs(pdist2([0,0,0],[c,p,gp])-pdist2([0,0,0],[c_new,p_new,gp_new])) < thresh  && ...
                    0 <= c_new && c_new <= 1 && 0 <= p_new && p_new <= 1 && 0<=gp_new && gp_new<=1 && c_new <= p_new && p_new<=gp_new  
                    r(i) = c_new;
                    s(i) = p_new;
                    t(i) = gp_new;
                    j = maxit;
                    continue;
              else 
                    c= c_new;
                    p = p_new;
                    gp = gp_new; 
              end
        end
    end 
end 
f_subs = [r;s;t];
end
 
function [ f_subs ] = Alt2D_Proj_II_graph(f_vec,maxit,f_3dproj)

N = length(f_vec);
n = N/3 ;
c_vec = f_vec(1:n);
p_vec = f_vec(n+1:2*n);
gp_vec = f_vec(2*n+1:3*n);

r = zeros(n,1);
s = zeros(n,1);
t = zeros(n,1);

for i = 1:n
    c = c_vec(i);
    p = p_vec(i);
    gp = gp_vec(i);

j =1; 

figure(1)
corners = [0 0 0;0 1 1;0 0 1;1 1 1]; 
scatter3(corners(:,1),corners(:,2),corners(:,3),'.');hold on;
plot3(corners(1:2,1),corners(1:2,2),corners(1:2,3),LineWidth=1, Color="red"); hold on;
plot3(corners(2:3,1),corners(2:3,2),corners(2:3,3),LineWidth=1, Color="red"); hold on;
plot3(corners(3:4,1),corners(3:4,2),corners(3:4,3),LineWidth=1, Color="red"); hold on;
plot3([1 0],[1 0],[1 0],LineWidth=1, Color="red"); hold on;
plot3([0 0],[0 0],[1 0],LineWidth=1, Color="red"); hold on;
plot3([0 1],[1 1],[1 1], LineWidth=1, Color="red"); hold on;
str = ['+','o','*','s','d','^','p','h'];
for ii=1:1
plot3(f_vec(1),f_vec(2),f_vec(3),str(ii),color="black",LineWidth=2); hold on 
plot3(f_3dproj(1),f_3dproj(2),f_3dproj(3),str(ii),color="red",LineWidth=2);hold on
end
xlabel('$C$','FontSize',30,Interpreter="latex");
ylabel('$P$','FontSize',30, Interpreter="latex");
zlabel('$GP$','FontSize',30, Interpreter="latex");
xlim([-0.5 1.5])
ylim([-0.5 1.5])
zlim([-0.5 1.5])

while (j < maxit)
    if j==1
    zppos = sort([0,gp,1]);
    gp = zppos(2);
    end
    %----------------------------------------------------------------------
    [c_new, p_new] = Projec_2D(c, p, 1);
    [c_new, gp_new] = Projec_2D(c_new, gp, 1); 
    [p_new, gp_new] = Projec_2D(p_new, gp_new, 1);
    j = j+1;
    %----------------------------------------------------------------------
    plot3(c_new,p_new,gp_new,str(ii),LineWidth=2);hold on
    pause(2); 
    if abs(pdist2([0,0,0],[c,p,gp])-pdist2([0,0,0],[c_new,p_new,gp_new])) < thresh
       r(i) = c_new;
       s(i) = p_new;
       t(i) = gp_new;
       j = maxit;
       continue;
    else 
       c= c_new;
       p = p_new;
       gp = gp_new; 
    end
end
hold off
end
f_subs = [r;s;t];
end

function [ f_subs ] = Alt2D_Proj_I(f_vec,maxit)


N = length(f_vec);
n = N/3; 

c_vec = f_vec(1:n);
p_vec = f_vec(n+1:2*n);
gp_vec = f_vec(2*n+1:3*n);

r = zeros(n,1);
s = zeros(n,1);
t = zeros(n,1);

for i = 1:n
    c = c_vec(i);
    p = p_vec(i);
    gp = gp_vec(i);

j =1; 
while (j < maxit)
    % Check if gp_parent-child observation already satisfies constraints.
    % If yes dont do projections and go to the next value 
    if 0 <= c && c <= 1 && 0 <= p && p <= 1 && 0<=gp && gp<=1 && c <= p && p<=gp   %Region inside the triangle
               r(i) = c;
               s(i) = p;
               t(i) = gp;
               j = maxit;
              continue;
    end  
    
    % Initialize gp, enforce it to be on 0<=gp<=1 
    if j==1
    zppos = sort([0,gp,1]);
    gp = zppos(2);
    end 
    %----------------------------------------------------------------------
    %     Fix gp,  project c and p
    %----------------------------------------------------------------------
     
    [c_new, p_new] = Projec_2D(c, p, gp);
    
    % Update c and p, then fix them to project onto 
    c = c_new;
    p = p_new; 

    [c_new, gp_new] = Projec_2D(c, gp, 1);

    % Update c and p, then fix them to project onto 
    c = c_new;
    gp = gp_new; 

    [p_new, gp_new] = Projec_2D(p, gp, 1);
    
    p =p_new;
    gp =gp_new;

    
    j = j+1;
end
end
f_subs = [r;s;t];
end
