clc
clear 
close all 

% This sript detect the SV on GP-> P-> C  Diploid case
% =====================================================================    
% =    True SV vectors and obsevration vectors Simulations    =
% =====================================================================
% this piece of code is to generate the vector f for GP, f for P and f for C

n = 10^5;        % length of the each signal (f vector)

%---- k_gp1 >= k_gp2
k_gp1 = round(0.05*n);          % number of SV in GP1, 5% of the signal
k_gp2 = round(0.05*n);             % number of SV in GP2, 5% of the signal

%k_p = round(0.5*k_gp);    % number of SV in P
%k_c = round(0.25*k_gp);   % number of SV in C

p = 0.7;              % percent of II copy SVs in p1 and p2, if p1 = 1, then the parent is completely homozygous, two copies 
%p2 = 0.0;            % percent of II copy SVs in p2
per = 0.9;            % similarity of two GPs If per = 1, gp1 and gp2 are identical
per2= 0.9;            % similarity percent of two ps 

%============================================================================================%
%------------------------- gp1 signal: Percentage from biology -------------------------------%
%============================================================================================%
numHomo_gp1 = floor(k_gp1*p);           % number of homozygous variants gp should have
numHetero_gp1 = k_gp1 - numHomo_gp1;    % number of heterozygous variants gp should have

numHomo_gp2 = floor(k_gp2*p);           % number of homozygous variants gp should have
numHetero_gp2 = k_gp2 - numHomo_gp2;    % number of heterozygous variants gp should have
%%% this information is to create another gp with some similarity 

numShare_Homo_gp   = min(numHomo_gp2,floor(numHomo_gp1*per));
numNoShare_Homo_gp = numHomo_gp2-numShare_Homo_gp;         % (1-persim)*numHomo
numShare_Hetero_gp = min(numHetero_gp2,floor(numHetero_gp1*per));         % number of heterozygous variants the parents should share  
numNoShare_Hetero_gp = numHetero_gp2-numShare_Hetero_gp;   % (1-persim)*numHetero

 
%------------------------------------------------------------------------%
q = randperm(n);  %1 to n
%================== y_gp1 and y_gp2 ====================
% Generate the parents heterozygous sparse vectors y_gp and y_p
y_gp1_true = zeros(n,1);
y_gp1_true(q(1:numHetero_gp1)) = 1;

%----- y_gp2 ----
y_gp2_true = zeros(n,1);
p2q  = randperm(floor(numShare_Hetero_gp));
y_gp2_true(q(p2q))         = 1;
Hetero_gp2 = numHetero_gp1 + numNoShare_Hetero_gp;
y_gp2_true(q((numHetero_gp1+1): Hetero_gp2)) = 1;

%----------- z_gp1 and z_gp2
z_gp1_true = zeros(n,1);
z_gp2_true = zeros(n,1);
p2q_New  = randperm(floor(numShare_Homo_gp))+ Hetero_gp2;  %this randomly permutes the indices that the parents will share

z_gp1_true(q(Hetero_gp2+1:Hetero_gp2 + numHomo_gp1))  = 1; % assigns the correct number of homozygous instances to parent f

z_gp2_true(q(p2q_New)) = 1;              % assigns the correct number of shared sites to parent m
HomoSig_end = Hetero_gp2 + numHomo_gp1 + numNoShare_Homo_gp;
z_gp2_true(q((Hetero_gp2 + numHomo_gp1 + 1): HomoSig_end)) = 1;


%===============================================================================================%
%-------------------------------- p signal: based on gp signal ---------------------------------%
%===============================================================================================%
[z_p1_true,y_p1_true]=f_true_simu(z_gp1_true,y_gp1_true,z_gp2_true,y_gp2_true); 

fprintf('=============== GP1 =============== \n')
fprintf('Number of homozygous SV in GP1 = %4.2f \n',sum(z_gp1_true))
fprintf('Number of heterozygous SV in GP1 = %4.2f \n',sum(y_gp1_true))
fprintf('Number of SVs in GP1 = %4.2f \n',sum(y_gp1_true) + sum(z_gp1_true))

fprintf('=============== GP2 =============== \n')
fprintf('Number of homozygous SV in GP2 = %4.2f \n',sum(z_gp2_true) )
fprintf('Number of heterozygous SV in GP2 = %4.2f \n',sum(y_gp2_true) )
fprintf('Number of SVs in GP2 = %4.2f \n',sum(y_gp2_true) + sum(z_gp2_true) )

fprintf('=============== P1 =============== \n')
fprintf('Number of homozygous SV in P1 = %4.2f \n',sum(z_p1_true) )
fprintf('Number of heterozygous SV in P1 = %4.2f \n',sum(y_p1_true) )
fprintf('Number of SVs in P1 = %4.2f \n', sum(y_p1_true) + sum(z_p1_true) )

%  for i=1:n
%     % 1-- if gp is homozygous, p must be either homozygous or hetro,
%     if z_gp_true(i)==1 
%            temp=randi(2);      % to more enforce homozygosity 
%            if temp ==1
%               y_p_true(i)=1;
%               z_p_true(i)=0;
%            else 
%               z_p_true(i)=1;  
%               y_p_true(i)=0;  
%            end 
%     end 
% 
%     % 2- P cannot have II copies of SV if GP dose not have at least one copy
%     %If both are 0 parent should have 0.
%     if z_gp_true(i) ==0 && y_gp_true(i)==0 
%         temp=randi(700);         % to more enforce sparsity  
%         if temp==1
%             z_p_true(i)=0;
%             y_p_true(i)=1;
%         else
%             z_p_true(i)=0;
%             y_p_true(i)=0;   
%         end 
%     end 
% 
%     %3- GP is hetero, P is either homo, hetro or non 
%     if y_gp_true(i)==1 
%       temp=randi(300);           % to more enforce sparsity  
%         if temp ==1
%             z_p_true(i)=1;
%             y_p_true(i)=0;
%         elseif temp ==2
%             y_p_true(i)=1;
%             z_p_true(i)=0;
%         else 
%             y_p_true(i)=0;
%             z_p_true(i)=0;
%         end
%     end 
% end 


%===============================================================================================%
%-------------------------------- c signal: based on p signal ----------------------------------%
%===============================================================================================%
%%% assuming p1>=p2 
numHomo_p1 = sum(z_p1_true);
numHetero_p1 = sum(y_p1_true);
p2= numHomo_p1/(numHetero_p1+numHomo_p1);     % percent of homo. SV in p2 
per2=0.9;                                     % similarity percent of two ps 

%--------- start simulating p2
k_p2 = sum(y_p1_true) + sum(z_p1_true);  %round(0.025*n);               % number of SV in P2, 2.5% of the signal 
numHomo_p2 = sum(z_p1_true);% floor(k_p2*p2);         % number of homozygous variants in p2 
numHetero_p2 = k_p2 - numHomo_p2;   % number of heterozygous variants p2 

%%% this information is to create another gp1 with some similarity 
numShare_Homo_p   = min(numHomo_p2,floor(numHomo_gp1*per2));
numNoShare_Homo_p = numHomo_p2-numShare_Homo_p;         % (1-persim)*numHomo
numShare_Hetero_p = min(numHetero_p2,floor(numHetero_p1*per2));         % number of heterozygous variants the parents should share  
numNoShare_Hetero_p = numHetero_p2-numShare_Hetero_p;   % (1-persim)*numHetero

%------------------------------------------------------------------------%
sv_y_p1= find(y_p1_true);
sv_y_p1 = sv_y_p1(randperm(length(sv_y_p1)));
sv_z_p1= find(z_p1_true);
sv_z_p1 = sv_z_p1(randperm(length(sv_z_p1)));

%==================  y_p2 ====================
% Generate the parents heterozygous sparse vectors y_gp and y_p
y_p2_true = zeros(n,1);
p2q  = randperm(floor(numShare_Hetero_p));
y_p2_true(sv_y_p1(p2q))= 1;

q2 = randperm(n);          %1 to n
q2 = setdiff(q2,sv_y_p1);
y_p2_true(q2(1: numNoShare_Hetero_p)) = 1;
%----------- z_p2
z_p2_true = zeros(n,1);
p2q_New  = randperm(floor(numShare_Homo_p));  
z_p2_true(sv_z_p1(p2q_New)) = 1;                % assigns the correct number of shared sites to parent m
q3 = randperm(n);          %1 to n
q3 = setdiff(q3,sv_z_p1);
z_p2_true(q3(1:numNoShare_Homo_p)) = 1;

fprintf('=============== P2 =============== \n')
fprintf('Number of homozygous SV in P2 = %4.2f \n',sum(z_p2_true) )
fprintf('Number of heterozygous SV in P2 = %4.2f \n',sum(y_p2_true) )
fprintf('Number of SVs in P2 = %4.2f \n', sum(y_p2_true) + sum(z_p2_true) )

[z_c_true,y_c_true]=f_true_simu(z_p1_true,y_p1_true,z_p2_true,y_p2_true); 

fprintf('=============== C =============== \n')
fprintf('Number of homozygous SV in C = %4.2f \n',sum(z_c_true) )
fprintf('Number of heterozygous SV in C = %4.2f \n',sum(y_c_true) )
fprintf('Number of SVs in P2 = %4.2f \n', sum(y_c_true) + sum(z_c_true) )



%%

%==============================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Observations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create observations based on y_c ~ Poisson ( lambda_c - eps) I*fc + epsi))
% Generate observations and A_c, A_p matrices
%==============================================================================
% Initialize lambda and epsilon values
lambda_c = 6;
lambda_p = 4;
lambda_gp =2;
epsilon = 0.3;
tau=1;

A_c_z = (2*lambda_c - epsilon)*speye(n);
A_c_y = (lambda_c - epsilon)*speye(n);
%A_c = [A_c_z, A_c_y];

A_p_z=(2*lambda_p - epsilon)*speye(n);
A_p_y=(lambda_p - epsilon)*speye(n);
%A_f = [A_f_z, A_f_y];

A_gp_z = (2*lambda_gp - epsilon)*speye(n);
A_gp_y = (lambda_gp - epsilon)*speye(n);
%A_m = [A_m_z, A_m_y];

A=blkdiag([A_c_z,A_c_y], [A_p_z,A_p_y], [A_gp_z,A_gp_y]);


% z_c_obs = poissrnd(A_c_z * z_c_true + epsilon * ones(n,1));
% z_p_obs = poissrnd(A_p_z * z_p_true + epsilon * ones(n,1));
% z_gp_obs = poissrnd(A_gp_z * z_gp_true + epsilon * ones(n,1));
% 
% y_c_obs = poissrnd(A_c_y * y_c_true + epsilon * ones(n,1));
% y_p_obs = poissrnd(A_p_y * y_p_true + epsilon * ones(n,1));
% y_gp_obs = poissrnd(A_gp_y * y_gp_true + epsilon * ones(n,1));

s_c = poissrnd(  A_c_z * z_c_true + A_c_y * y_c_true + epsilon * ones(n,1));
s_p = poissrnd(  A_p_z * z_p1_true + A_p_y * y_p1_true + epsilon * ones(n,1));
s_gp = poissrnd(  A_gp_z * z_gp1_true + A_gp_y * y_gp1_true + epsilon * ones(n,1));

f = [z_c_true; y_c_true; z_p1_true; y_p1_true; z_gp1_true; y_gp1_true];        % True Signal
s = [s_c;s_p;s_gp];     % Observed Data

%%
N = length(f);
n = N/6;

% Setup function handles for computing A and A^T:
AT  = @(x) A'*x;
Ax   = @(x) A*x;

maxiter = 1000;
tolerance = 1e-8;
% verbose is how often it prints to the screen
verbose = 500;

% Simple initialization:
% AT(y) rescaled to a least-squares fit to the mean intensity

finit = (sum(sum(s)).*numel(AT(s)))...
    ./(sum(sum(AT(s))) .*sum(sum((AT(ones(size(s)))))))...
    .*AT(s);


[fhatSPIRAL_2p1c_dip, iterationsSPIRAL, objectiveSPIRAL,...
    reconerrorSPIRAL, cputimeSPIRAL] ...
    = SPIRALTAP_Dip_2P1C(s,Ax,tau,...
    'maxiter',maxiter,...
    'Initialization',finit,...
    'AT',AT,...
    'miniter',5,...
    'stopcriterion',3,...
    'tolerance',tolerance,...
    'alphainit',1,...
    'alphamin', 1e-30,...
    'alphamax', 1e30,...
    'alphaaccept',1e30,...
    'logepsilon',1e-10,...
    'saveobjective',1,...
    'savereconerror',1,...
    'savecputime',1,...
    'savesolutionpath',0,...
    'truth',f,...
    'verbose',verbose);


fhatSPIRAL_dip_c_z = fhatSPIRAL_2p1c_dip(1:n);
fhatSPIRAL_dip_c_y = fhatSPIRAL_2p1c_dip(n+1:2*n);
fhatSPIRAL_dip_p_z = fhatSPIRAL_2p1c_dip(2*n+1:3*n);
fhatSPIRAL_dip_p_y = fhatSPIRAL_2p1c_dip(3*n+1:4*n);
fhatSPIRAL_dip_gp_z = fhatSPIRAL_2p1c_dip(4*n+1:5*n);
fhatSPIRAL_dip_gp_y = fhatSPIRAL_2p1c_dip(5*n+1:6*n);


thresh = linspace(-.001,1.001,100); % make thresholding vector for ROC curve


%% 
%--------------------- 2p1C Diploid Method ROC -------------------------

f_true = f;
f_spiral = [fhatSPIRAL_dip_c_z;fhatSPIRAL_dip_c_y;fhatSPIRAL_dip_p_z;fhatSPIRAL_dip_p_y;fhatSPIRAL_dip_gp_z;fhatSPIRAL_dip_gp_y];

figure(1)
[FPRv_cpgp,TPRv_cpgp,PPV_cp] = ROC_CURVE(f_spiral, f_true, thresh);
%ROC_data = roc_curve(f_spiral, f_true, 1);
AUC = 1-trapz(1-TPRv_cpgp,FPRv_cpgp);
[recall,accuracy,precision] = metrics(f_true,f_spiral) ; 

figure(2)  
plot(FPRv_cpgp,TPRv_cpgp,'blue','LineWidth',3); hold on;

plot([0,1],[0,1],'r-.','LineWidth',1); hold off;
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
title("ROC Curve - AUC = "+ AUC)

%%
%[AUC] = perfcurve(f_spiral>= 0.5,f_true,'virginica')
function ROC_data = roc_curve(class_1, class_2, dispp, dispt)
    % Setting default parameters and detecting errors
    if(nargin<4), dispt = 1;    end
    if(nargin<3), dispp = 1;    end
    if(nargin<2), error('Params "class_1" or "class_2" are not indicated.'); end
    class_1 = class_1(:);
    class_2 = class_2(:);
    
    % Calculating the threshold values between the data points
    s_data = unique(sort([class_1; class_2]));          % Sorted data points
    s_data(isnan(s_data)) = [];                 % Delete NaN values
    d_data = diff(s_data);                      % Difference between consecutive points
    if(isempty(d_data)), error('Both class data are the same!'); end
    d_data(length(d_data)+1,1) = d_data(length(d_data));% Last point
    thres(1,1) = s_data(1) - d_data(1);                 % First point
    thres(2:length(s_data)+1,1) = s_data + d_data./2;   % Threshold values
        
    % Calculating the sensibility and specificity of each threshold
    curve = zeros(size(thres,1),2);
    distance = zeros(size(thres,1),1);
    for id_t = 1:1:length(thres)
        TP = length(find(class_2 >= thres(id_t)));    % True positives
        FP = length(find(class_1 >= thres(id_t)));    % False positives
        FN = length(find(class_2 < thres(id_t)));     % False negatives
        TN = length(find(class_1 < thres(id_t)));     % True negatives
        
        curve(id_t,1) = TP/(TP + FN);   % Sensitivity
        curve(id_t,2) = TN/(TN + FP);	% Specificity
        
        % Distance between each point and the optimum point (0,1)
        distance(id_t)= sqrt((1-curve(id_t,1))^2+(curve(id_t,2)-1)^2);
    end
    
    % Optimum threshold and parameters
    [~, opt] = min(distance);
    TP = length(find(class_2 >= thres(opt)));    % No. true positives
    FP = length(find(class_1 >= thres(opt)));    % No. false positives 
    FN = length(find(class_2 < thres(opt)));     % No. false negatives                                 
    TN = length(find(class_1 < thres(opt)));     % No. true negatives       
    
    % Output parameters
    param.Threshold = thres(opt);               % Optimum threshold position
    param.Sensi = curve(opt,1);                 % Sensitivity
    param.Speci = curve(opt,2);                 % Specificity
    param.AROC  = abs(trapz(1-curve(:,2),curve(:,1))); % Area under curve
    param.Accuracy = (TP+TN)/(TP+TN+FP+FN);     % Aaccuracy
    param.PPV   = TP/(TP+FP);                   % Positive predictive value
    param.NPV   = TN/(TN+FN);                   % Negative predictive value
    param.FNR   = FN/(FN+TP);                   % False negative rate
    param.FPR   = FP/(FP+TN);                   % False positive rate
    param.FDR   = FP/(FP+TP);                   % False discovery rate
    param.FOR   = FN/(FN+TN);                   % False omission rate
    param.F1_score = 2*TP/(2*TP+FP+FN);         % F1 score
    param.MCC   = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));  % Matthews correlation coefficient
    param.BM    = param.Sensi+param.Speci-1;    % Informedness
    param.MK    = param.PPV+param.NPV-1;        % Markedness
    
    param.TP = TP;    % No. true positives
    param.FP = FP;    % No. false positives 
    param.FN = FN;    % No. false negatives                                 
    param.TN = TN;    % No. true negatives  
    
    % Plotting if required
    if(dispp == 1)
        fill_color = [11/255, 208/255, 217/255];
        fill([1-curve(:,2); 1], [curve(:,1); 0], fill_color,'FaceAlpha',0.5);
        hold on; plot(1-curve(:,2), curve(:,1), '-b', 'LineWidth', 2);
        hold on; plot(1-curve(opt,2), curve(opt,1), 'or', 'MarkerSize', 10);
        hold on; plot(1-curve(opt,2), curve(opt,1), 'xr', 'MarkerSize', 12);
        hold off; axis square; grid on; xlabel('1 - specificity'); ylabel('sensibility');
        title(['AROC = ' num2str(param.AROC)]);
    end
    
    % AROC warning
    if param.AROC < 0.5
        warning('Since AROC is less than 0.5, you should swap the classes: roc_curve(class_2,class_1).');
    end
    
    % Log screen parameters if required
    if(dispt == 1)
        fprintf('\n ROC CURVE PARAMETERS\n');
        fprintf(' ------------------------------\n');
        fprintf('  - Distance:     %.4f\n', distance(opt));
        fprintf('  - Threshold:    %.4f\n', param.Threshold);
        fprintf('  - Sensitivity:  %.4f\n', param.Sensi);
        fprintf('  - Specificity:  %.4f\n', param.Speci);
        fprintf('  - AROC:         %.4f\n', param.AROC);
        fprintf('  - Accuracy:     %.4f\n', param.Accuracy);
        fprintf('  - PPV:          %.4f\n', param.PPV);
        fprintf('  - NPV:          %.4f\n', param.NPV);
        fprintf('  - FNR:          %.4f\n', param.FNR);
        fprintf('  - FPR:          %.4f\n', param.FPR);
        fprintf('  - FDR:          %.4f\n', param.FDR);
        fprintf('  - FOR:          %.4f\n', param.FOR);
        fprintf('  - F1 score:     %.4f\n', param.F1_score);
        fprintf('  - MCC:          %.4f\n', param.MCC);
        fprintf('  - BM:           %.4f\n', param.BM);
        fprintf('  - MK:           %.4f\n', param.MK);
        fprintf(' \n');
    end
    
    % Assinging parameters and curve data
    ROC_data.param = param;
    ROC_data.curve = curve;
end 

function [ FPRv, TPRv, PPV ] = ROC_CURVE( fhat_recon, f_true, thresh)

% NOVroc_gen_PPV takes in SPIRAL reconstruction, true signal, and thresholding
% values and outputs the false and true positive rates to build ROC and PPV curves.

% Initialize vectors to store True positives and False positives values
T_pv = zeros(length(thresh),1);
F_pv = zeros(length(thresh),1);

T_nv = zeros(length(thresh),1);
F_nv = zeros(length(thresh),1);

TPRv = zeros(length(thresh),1);
FPRv = zeros(length(thresh),1);

PPV = zeros(length(thresh),1);

% Determine size of true and reconstructed signals
n = length(fhat_recon);

for i = 1:length(thresh)
    
    f_thresh = fhat_recon >= thresh(i);
   
    
    T_p = 0;
    F_p = 0;
    
    T_n = 0;
    F_n = 0;   
    
    for j = 1:n
        
    if f_true(j)==1 && f_thresh(j)==1 %&& fhatSPIRAL_uc(j)~=0
        T_p = T_p + 1;
        
    elseif f_true(j)==0 && f_thresh(j)==1 %&&fhatSPIRAL_uc(j)~=0
        F_p = F_p + 1;
        
    elseif  f_true(j)==0 && f_thresh(j)==0 %&& fhatSPIRAL_uc(j)==0
        T_n = T_n + 1;   
        
    elseif f_true(j)==1 && f_thresh(j)==0 %&& fhatSPIRAL_uc(j)==0 
        F_n = F_n + 1;    
    end
    
    end
    
    T_pv(i) = T_p;
    F_pv(i) = F_p;
    
    T_nv(i) = T_n;
    F_nv(i) = F_n;
    
    
    TPRv(i) = T_p/(T_p + F_n);
    FPRv(i) = F_p/(F_p + T_n);
    
    
    %Positive Predictive Value: (TP)/(TP+FP)
    PPV(i) = T_p/(T_p+F_p);
    
end

end

function [Recall,Accur,precision] = metrics(f_true,fhat_recon)
thresh=0.3; 
f_thresh = fhat_recon >= thresh;
TP = sum((f_thresh == 1) & (f_true == 1)); 
FP = sum((f_thresh == 1) & (f_true == 0));
TN = sum((f_thresh == 0) & (f_true == 0));
FN = sum((f_thresh == 0) & (f_true == 1));
Accur = (TP+TN)./(TP+FP+TN+FN); 
Recall = TP/(TP + FN);
precision = TP/(TP + FP);
end

function [z_p1_true,y_p1_true]=f_true_simu(z_gp1_true,y_gp1_true,z_gp2_true,y_gp2_true)
y_p1_true = zeros(length(z_gp2_true),1);
z_p1_true = zeros(length(z_gp2_true),1);
n = length(z_gp2_true); 
for i=1:n 
   if z_gp1_true(i)==1 && y_gp1_true(i)==0 && z_gp2_true(i)==1  && y_gp2_true(i)==0  
           z_p1_true(i)=1;
           y_p1_true(i)=0; 
   elseif z_gp1_true(i)==1 && y_gp1_true(i)==0 && z_gp2_true(i)==0  && y_gp2_true(i)==1
              temp=randi(2);           % To more enforce homozygosity 
           if temp ==1
              z_p1_true(i)=0;
              y_p1_true(i)=1;
           else
              z_p1_true(i)=1;
              y_p1_true(i)=0;
           end
   elseif z_gp1_true(i)==1 && y_gp1_true(i)==0 && z_gp2_true(i)==0  && y_gp2_true(i)==0
              z_p1_true(i)=0;
              y_p1_true(i)=1;
   elseif z_gp1_true(i)==0 && y_gp1_true(i)==1  && z_gp2_true(i)==1 && y_gp2_true(i)==0
              temp=randi(2);           
           if temp ==1
              z_p1_true(i)=0;
              y_p1_true(i)=1;
           else
              z_p1_true(i)=1;
              y_p1_true(i)=0;
           end
   elseif  z_gp1_true(i)==0 && y_gp1_true(i)==1 && z_gp2_true(i)==0 && y_gp2_true(i)==1
              temp=randi(4); 
           if temp ==4
              z_p1_true(i)=1;
              y_p1_true(i)=0;
           elseif temp ==2 
              z_p1_true(i)=0;
              y_p1_true(i)=0;
           else 
              z_p1_true(i)=0;
              y_p1_true(i)=1;
           end 
   elseif  z_gp1_true(i)==0 && y_gp1_true(i)==1 && z_gp2_true(i)==0  && y_gp2_true(i) == 0      
              temp=randi(2);                % To more enforce homozygosity
           if temp ==1                               
              z_p1_true(i)=0;                
              y_p1_true(i)=0;                
           else                              
              z_p1_true(i)=0;
              y_p1_true(i)=1;
           end   
    elseif z_gp1_true(i)==0 && y_gp1_true(i)==0  && z_gp2_true(i)==1 && y_gp2_true(i)==0
              z_p1_true(i)=0;
              y_p1_true(i)=1;
   elseif  z_gp1_true(i)==0 && y_gp1_true(i)==0 && z_gp2_true(i)==0 && y_gp2_true(i)==1
              temp=randi(2);                % To more enforce homozygosity
           if temp ==1                               
              z_p1_true(i)=0;                
              y_p1_true(i)=0; 
           else                              
              z_p1_true(i)=0;
              y_p1_true(i)=1;
           end 
   else 
          z_p1_true(i)=0;
          y_p1_true(i)=0;
   end
end 
end



