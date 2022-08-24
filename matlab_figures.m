%% [1]  Precision-recall curve 
clc 
close all 
clear
%--------------Inputs----------
true_la = round(rand(100,1));    % true labels
pred1 = round(rand(100,1));      % predictions of model 1
pred2 = round(rand(100,1));      % predictions of model 2 
%------------------------------
figure(1)
[X1,Y1,T1,AUC1] = perfcurve(true_la,pred1,'1','XCrit', 'tpr', 'YCrit', 'prec');
[X2,Y2,T2,AUC2] = perfcurve(true_la,pred2,'1','XCrit', 'tpr', 'YCrit', 'prec');
Y1(1)=1;
Y2(1)=1;
Y1(end)=0;
Y2(end)=0;
Area_recall_3d=trapz(-Y1,X1);
Area_recall_2d=trapz(-Y2,X2); 
plot(X1,Y1,'blue','LineWidth',3) ;hold on;
plot(X2,Y2,'red','LineWidth',3,'LineStyle','--') ;hold on;
plot([0,1],[1,0],'r-.','LineWidth',1);    
xlabel('Recall',FontSize=15)
ylabel('Precision',FontSize=15)
dim1 = [.2 .3 .3 .1];
dim2 = [.2 .3 .3 .01];
str1 = strcat('AUC_{3D}= ',num2str(round(Area_recall_2d,3)));
str2 = strcat('AUC_{2D}= ',num2str(round(Area_recall_3d,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
legend("3D Projections","2D Projections","Location","southwest")
title('PR Curve for Reconstructed Child Signal')

%% [2]  ROC Curve  
clc
close all 
clear 
%--------------Inputs----------
true_la = round(rand(100,1));    % true labels
pred1 = round(rand(100,1));      % predictions of model 2
pred2 = round(rand(100,1));      % predictions of model 3 
%------------------------------
figure(2)
[X1,Y1,T1,AUC1] = perfcurve(true_la,pred1,'1');     % 3d
[X2,Y2,T2,AUC2] = perfcurve(true_la,pred2,'1');   % 2d 
AreaUnerCurve_3d =1+trapz(1-Y1,X1);
AreaUnerCurve_2d =1+trapz(1-Y2,X2);
plot(X1,Y1,'blue','LineWidth',3) ;hold on;
plot(X2,Y2,'red','LineWidth',3) ;hold on;
plot([0,1],[0,1],'r-.','LineWidth',1); hold off;   
xlabel('False positive rate') 
ylabel('True positive rate')
dim1 = [.2 .5 .3 .1];
dim2 = [.2 .5 .3 .01];
str1 = strcat('AUC_{3D}= ',num2str(round(AUC1,3)));
str2 = strcat('AUC_{2D}= ',num2str(round(AUC2,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
legend("3D Projections","2D Projections","Location","southwest")
title('Grandparent Signal')

%% [3]  Random points within a linear feasible reagion   
clc
close all 
clear 
n=10000;     
amp=2;
phase=0.5;
point = amp*rand(n,3)- phase;

figure(3)
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

%% [4]  Stem plots 

figure(4)
dist_from_2d = round(rand(300,1)); 
dist_from_3d = round(rand(300,1)); 

num_poi=300;
stem( dist_from_2d(1:num_poi), 'filled','red');    hold on;
stem( dist_from_3d(1:num_poi), 'filled','blue')    

xlabel('Point Number','FontSize',30,Interpreter="latex");
ylabel('Euclidean Distance to Original Point','FontSize',20, Interpreter="latex");
legend("2D","3D") 
ylim([-0.01 2.5]) 
hold off 

%% [5]  Distributions - two overlaping distributions  

%--------------Inputs----------
true_la = round(rand(100,1));    % true labels
pred1 = round(10*rand(100,1));      % predictions of model 1
pred2 = round(10*rand(100,1));      % predictions of model 2 

%-------------------
[H1,e1] =hist(pred1);
[H2,e2] = hist(pred2);
h1 = H1/sum(H1);
h2 = H2/sum(H2); 
%----------------
figure(5)
histogram(pred1,'Normalization','pdf','FaceAlpha',0.2); hold on
histogram(pred2,'Normalization','pdf','FaceAlpha',0.2); hold off 
legend("Zeros","Ones")
xlabel("Coverage")
title ("Predicted - 0.5 threasholding") 

%% [6] Confusion Matix   
%--------------Inputs----------
true_la = round(rand(100,1));    % true labels
pred1 = round(rand(100,1));      % predictions of model 1
pred2 = round(rand(100,1));      % predictions of model 2 

figure(6)
cm = confusionchart(double(pred1),double(pred2))
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';








