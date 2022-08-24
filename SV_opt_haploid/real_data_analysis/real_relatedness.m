clc
clear 
close all 

True_sig = readtable("real_data/5X/17CEU_5X_SPIRAL_trueSignal.csv"); 
Observs = readtable("real_data/5X/17CEU_5X_ObsSignal.csv");  
% 1:NA12881: GM Paternal Grandmother               
% 2:NA12881: C daughter
% 3:NA12886: C Son
% 4:NA12883: C Son
% 5:NA12879: C daughter
% 6:NA12888: C Son
% 7:NA12892: GM Maternal Grandmother
% 8:NA12889: GP Paternal Grandfather  <---------------
% 9:NA12885: C daughter 
% 10:NA12882: C Son 
% 11:NA12880: C daughter 
% 12:NA12877: P Father            <---------------
% 13:NA12891: GP Maternal Grandfather  
% 14:NA12878: P Mother 
% 15:NA12884: C Son
% 16:NA12893: C Son 
% 17:NA12887: C daughter
%%
f_c = True_sig{:,11};
f_c(f_c>0)=1;
f_p = True_sig{:,12};
f_p(f_p>0)=1;
f_gp = True_sig{:,8};
f_gp(f_gp>0)=1;

y_c = Observs{:,11};
y_p = Observs{:,12};
y_gp = Observs{:,8};


%%
thres = 7;
ind_c_0 =find(y_c ==0); 
ind_c_20 =find(y_c >thres); 
ind_c = [ind_c_0;ind_c_20];
y_c = y_c(ind_c); 
y_p = y_p(ind_c); 
y_gp = y_gp(ind_c); 
f_c = f_c(ind_c);
f_p = f_p(ind_c);
f_gp = f_gp(ind_c);

ind_p_0 =find(y_p ==0); 
ind_p_20 =find(y_p >thres); 
ind_p = [ind_p_0;ind_p_20];
y_c = y_c(ind_p); 
y_p = y_p(ind_p); 
y_gp = y_gp(ind_p); 
f_c = f_c(ind_p);
f_p = f_p(ind_p);
f_gp = f_gp(ind_p);

ind_gp_0 =find(y_gp ==0); 
ind_gp_20 =find(y_gp >thres); 
ind_gp = [ind_gp_0;ind_gp_20];
y_c = y_c(ind_gp); 
y_p = y_p(ind_gp); 
y_gp = y_gp(ind_gp); 
f_c = f_c(ind_gp);
f_p = f_p(ind_gp);
f_gp = f_gp(ind_gp);



length(y_c)
length(f_p)

length(find(f_c==1))
length(find(f_p==1))
length(find(f_gp==1))

length(find(y_c > 20))
length(find(y_p > 20))
length(find(y_gp > 20))



%%
close all 
figure(1)

subplot(1,3,1)
histogram(y_gp)
xlim([1 8])
ylim([1 70000])

title("Histogram of $y_{gp}$",Interpreter="latex",FontSize=20)
xlabel("Number of Fregments")
ylabel("Counts")
subplot(1,3,2)
histogram(y_p)
xlim([1 8])
ylim([1 70000])

title("Histogram of $y_{p}$",Interpreter="latex",FontSize=20)
xlabel("Number of Fregments")
ylabel("Counts")

subplot(1,3,3)
histogram(y_c)
xlim([1 8])
ylim([1 70000])

title("Histogram of $y_{c}$",Interpreter="latex",FontSize=20)
xlabel("Number of Fregments")
ylabel("Counts")



%%
n = length(y_gp);

% define parameters 
lambda_gp = 5;
lambda_p = 5;
lambda_c = 5;
epsilon = 0.2;
tau = 10.5; 
beta = 1; 
reg= [tau*beta^2,tau*beta,tau]; 
% define the coverage matrix A and calculate the
% observation y for each individual

A_gp = (lambda_gp - epsilon)*speye(n);
A_p = (lambda_p - epsilon)*speye(n);
A_c = (lambda_c - epsilon)*speye(n); 

%%% end of the section %%%

% =========================================================================    
% ==========    Set up GP/P/C Novel Method reconstruction  ================
% =========================================================================

 f = [f_c; f_p; f_gp];    % concatinate true signals
 y = [y_c; y_p; y_gp];    % concatinate observations

 N = length(f);
 n = N/3;

% set up the block diagonal matrix A
A = sparse(3*n, 3*n); 
A(1:n,1:n) = A_c;
A(n+1:2*n,n+1:2*n) = A_p;
A(2*n+1:3*n, 2*n+1:3*n) = A_gp;

% Setup function handles for computing A and A^T:
AT  = @(x) A'*x;
A   = @(x) A*x;

%%%%%%%%%%%%%%%% SPIRAL-TAP Paremeters %%%%%%%%%%%%%%%%
% set maximum number of iterations, tol, and when to print to screen
maxiter = 90;
tolerance = 1e-8;
verbose = 100;

%------------- initialize the f vector ------------------------------
        % Simple initialization:
        % AT(y) rescaled to a least-squares fit to the mean intensity
finit = (sum(sum(y)).*numel(AT(y)))...
        ./(sum(sum(AT(y))) .*sum(sum((AT(ones(size(y))))))).*AT(y);

        [fhatSPIRAL_1gp1p1c, iterationsSPIRAL, objectiveSPIRAL,...
            reconerrorSPIRAL, cputimeSPIRAL] ...
            = SPIRALTAP_1GP1P1C_Hap_3d(y,A,reg,...
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
            'logepsilon',epsilon,...
            'saveobjective',1,...
            'savereconerror',1,...
            'savecputime',1,...
            'savesolutionpath',0,...
            'truth',f,...
            'verbose',verbose);


% separate reconstruction for each true signal   3d 
  fhatSPIRAL_c = fhatSPIRAL_1gp1p1c(1:n);
  fhatSPIRAL_p = fhatSPIRAL_1gp1p1c(n+1:2*n);
  fhatSPIRAL_gp = fhatSPIRAL_1gp1p1c(2*n+1:3*n);
  
  %---------------------------------------------------------
  % THRESHOLDING OF SIGNALS
  %---------------------------------------------------------
thresh = linspace(-.001,1.000,100); % make thresholding vector for ROC curve

%  Plot Results of when three signals together      
f_true = [f_c; f_p; f_gp];  
f_spiral = [fhatSPIRAL_c;fhatSPIRAL_p;fhatSPIRAL_gp];  

%%%%%%%%% Recall-precision Plot %%%%%%%%%%%%%%%
%%
close all 
figure(1)

[X1,Y1,T1,AUC1] = perfcurve(f_true,f_spiral,'1','XCrit', 'tpr', 'YCrit', 'prec');
Y1(1)=1;
Y2(1)=1;
Y1(end)=0;
Y2(end)=0;

Area_recall_3d=trapz(-Y1,X1);

plot(X1,Y1,'blue','LineWidth',3) ;hold on;
plot([0,1],[1,0],'r-.','LineWidth',1); hold off;   
xlabel('Recall')
ylabel('Precision')
dim1 = [.2 .5 .3 .1];

str1 = strcat('AUC_{3D}= ',num2str(round(Area_recall_3d,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
%legend("3D Projections","2D Projections","Location","best")
title('All Signals')

%=========================================

f_c_true = [f_c]; %; f_p; f_gp];  
f_c_spiral = [fhatSPIRAL_c];%fhatSPIRAL_p;fhatSPIRAL_gp];  

figure(2)
subplot(1,3,1);

[X1,Y1,T1,AUC1] = perfcurve(f_c_true,f_c_spiral,'1','XCrit', 'tpr', 'YCrit', 'prec');
Y1(1)=1;
Y2(1)=1;
Y1(end)=0;
Y2(end)=0;

Area_recall_3d=trapz(-Y1,X1);
plot(X1,Y1,'blue','LineWidth',3) ;hold on;
plot([0,1],[1,0],'r-.','LineWidth',1); hold off;   
xlabel('Recall')
ylabel('Precision')
dim1 = [.24 .3 .3 .1];
str1 = strcat('AUC= ',num2str(round(Area_recall_3d,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
%legend("3D Projections","2D Projections","Location","best")
title('PR Curve for Reconstructed Child Signal')

%===============================
f_p_true = [f_p]; %; f_p; f_gp];  
f_p_spiral = [fhatSPIRAL_p];%fhatSPIRAL_p;fhatSPIRAL_gp];  


subplot(1,3,2);
[X1,Y1,T1,AUC1] = perfcurve(f_p_true,f_p_spiral,'1','XCrit', 'tpr', 'YCrit', 'prec');
Y1(1)=1;
Y2(1)=1;
Y1(end)=0;
Y2(end)=0;
%Area_recall_3d=trapz(-Y1(2:end),X1(2:end));
%Area_recall_2d=trapz(-Y2(2:end),X2(2:end)); 
Area_recall_3d=trapz(1-Y1,X1);
plot(X1,Y1,'blue','LineWidth',3) ;hold on;
plot([0,1],[1,0],'r-.','LineWidth',1); hold off;   
xlabel('Recall')
ylabel('Precision')
dim1 = [.53 .3 .3 .1];
str1 = strcat('AUC= ',num2str(round(Area_recall_3d,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
%legend("3D Projections","2D Projections","Location","best")
title('PR Curve for Reconstructed Parent Signal')
%hold off 

%===============================
f_gp_true = [f_gp]; %; f_p; f_gp];  
f_gp_spiral = [fhatSPIRAL_gp];%fhatSPIRAL_p;fhatSPIRAL_gp];  


%figure(8)
subplot(1,3,3);
[X1,Y1,T1,AUC1] = perfcurve(f_gp_true,f_gp_spiral,'1','XCrit', 'tpr', 'YCrit', 'prec');
Y1(1)=1;
Y2(1)=1;
Y1(end)=0;
Y2(end)=0;
%Area_recall_3d=trapz(-Y1(2:end),X1(2:end));
%Area_recall_2d=trapz(-Y2(2:end),X2(2:end)); 
Area_recall_3d=trapz(1-Y1,X1);
plot(X1,Y1,'blue','LineWidth',3) ;hold on;
plot([0,1],[1,0],'r-.','LineWidth',1); hold off;   
xlabel('Recall')
ylabel('Precision')
dim1 = [.8 .3 .3 .1];
str1 = strcat('AUC= ',num2str(round(Area_recall_3d,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
%legend("3D Projections","2D Projections","Location","best")
title('PR Curve for Reconstructed Grandparent Signal')
hold off 

%%  ROC Curves 
close all
thresh = linspace(-.001,1.001,100); % make
[FPRv_cpgp,TPRv_cpgp,PPV_cp] = ROC_CURVE(f_spiral, f_true, thresh);

%ROC_data = roc_curve(f_spiral, f_true, 1);
AUC = 1-trapz(1-TPRv_cpgp,FPRv_cpgp);

%[recall,accuracy,precision] = metrics(f_true,f_spiral) ; 

figure(1)  
plot(FPRv_cpgp,TPRv_cpgp,'red','LineWidth',3); hold on;

plot([0,1],[0,1],'r-.','LineWidth',1); hold off;
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
dim1 = [.7 .3 .3 .1];
str1 = strcat('AUC= ',num2str(round(AUC,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');

%title("ROC Curve - AUC = "+ AUC)

%%%%%%%%%%%% Split figures %%%%%%%%%%%%%%%%

%%%%%%%%%%%% Child %%%%%%%%%%%%%%%%
figure(2)
f_c_true = [f_c]; %; f_p; f_gp];  
f_c_spiral = [fhatSPIRAL_c];%fhatSPIRAL_p;fhatSPIRAL_gp];  
thresh = linspace(-.001,1.001,100); % make
[FPRv_cpgp,TPRv_cpgp,PPV_cp] = ROC_CURVE(f_c_spiral, f_c_true, thresh);

%ROC_data = roc_curve(f_spiral, f_true, 1);
AUC = 1-trapz(1-TPRv_cpgp,FPRv_cpgp);

%[recall,accuracy,precision] = metrics(f_true,f_spiral) ; 

subplot(1,3,1)
plot(FPRv_cpgp,TPRv_cpgp,'red','LineWidth',3); hold on;

plot([0,1],[0,1],'r-.','LineWidth',1); hold off;
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
dim1 = [.21 .3 .3 .1];
str1 = strcat('AUC= ',num2str(round(AUC,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
title("Child Signal ")

%%%%%%%%%%%% Parent %%%%%%%%%%%%%%%%
f_p_true = [f_p]; %; f_p; f_gp];  
f_p_spiral = [fhatSPIRAL_p];%fhatSPIRAL_p;fhatSPIRAL_gp];  

thresh = linspace(-.001,1.001,100); % make
[FPRv_cpgp,TPRv_cpgp,PPV_cp] = ROC_CURVE(f_p_spiral, f_p_true, thresh);

%ROC_data = roc_curve(f_spiral, f_true, 1);
AUC = 1-trapz(1-TPRv_cpgp,FPRv_cpgp);

%[recall,accuracy,precision] = metrics(f_true,f_spiral) ; 

subplot(1,3,2)
plot(FPRv_cpgp,TPRv_cpgp,'red','LineWidth',3); hold on;

plot([0,1],[0,1],'r-.','LineWidth',1); hold off;
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
%legend("2D","3D")
dim1 = [.53 .3 .3 .1];
str2 = strcat('AUC= ',num2str(round(AUC,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
title("Parent Signal ")

%%%%%%%%%%%% Grand parent %%%%%%%%%%%%%%%%
f_gp_true = [f_gp]; %; f_p; f_gp];  
f_gp_spiral = [fhatSPIRAL_gp];%fhatSPIRAL_p;fhatSPIRAL_gp];  
thresh = linspace(-.001,1.001,100); % make
[FPRv_cpgp,TPRv_cpgp,PPV_cp] = ROC_CURVE(f_gp_spiral, f_gp_true, thresh);

%ROC_data = roc_curve(f_spiral, f_true, 1);
AUC = 1-trapz(1-TPRv_cpgp,FPRv_cpgp);

%[recall,accuracy,precision] = metrics(f_true,f_spiral) ; 

subplot(1,3,3)
plot(FPRv_cpgp,TPRv_cpgp,'red','LineWidth',3); hold on;

plot([0,1],[0,1],'r-.','LineWidth',1); hold off;
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
dim1 = [.8 .3 .3 .1];
str1 = strcat('AUC= ',num2str(round(AUC,3)));
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
title("Grandparent Signal ")


%% Helpful functions

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


%%%%%%%%%%%%%%%%%%%%%%%%%% SPIRAL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, varargout] = SPIRALTAP_1GP1P1C_Hap_3d(y, A, tau, varargin)
% ==== Set default/initial parameter values ====
% ---- All Methods ----
% =============================================================================
% =        SPIRAL:  Sparse Poisson Intensity Reconstruction Algorithms        =
% =                                Version 1.0                                =
% =============================================================================
% =    Copyright 2009, 2010                                                   =
% =    Zachary T. Harmany*, Roummel F. Marcia**, Rebecca M. Willett*          =
% =        *  Department of Electrical and Computer Engineering               =
% =           Duke University                                                 =
% =           Durham, NC 27708, USA                                           =
% =       **  School of Natural Sciences                                      =
% =           University of California, Merced                                =
% =           Merced, CA 95343, USA                                           =
% =                                                                           =
% =    Corresponding author: Zachary T. Harmany (zth@duke.edu)                =
% ============================================================================= 
%
% =============================================================================
% =                               Documentation                               =
% =============================================================================
% Syntax:
%   [x, optionalOutputs] = SPIRALTAP(y, A, tau, optionalInputs)
% 
%   More details and supporting publications are 
%   available on the SPIRAL Toolbox homepage
%   http://www.ee.duke.edu/~zth/spiral/
% 
% =============================================================================
% =                                  Inputs                                   =
% =============================================================================
% Required Inputs:
%	y               Degraded observations.  For this documenation, we say
%                   that y has m total elements.
%
%   A               Sensing / observation matrix.  For this documentation, 
%                   we say that A is an m x n matrix.  A could also be a
%                   function call A() that computes matrix-vector products
%                   such that A(x) = A*x where x has n total elements.  In 
%                   this case one must also specify a function call AT() 
%                   using the 'AT' option that computes matrix-vector 
%                   products with the adjoint of A such that AT(x) = A'*x.  
%
%   tau             Regularization parameter that trades off the data fit
%                   (negative log-likelihood) with the regularization.
%                   The regularization parameter can either be a
%                   nonnegative real scalar or (for all methods except the
%                   total variation penalty) have n nonnegative real 
%                   elements which allows for nonuniform penalization 
%                   schemes. 
%           
% Optional Inputs:
% If one were to only input y, A, and tau into the algorithm, there are 
% many necessary assumptions as to what to do with the inputs.  By default
% SPIRAL assumes that:
%   - y contains Poisson realizations of A*f, where f is the true underlying
%     signal (to be estimated),
%   - the penalty is the l_1 norm of x (i.e., we promote sparsity in the 
%     canonical basis.
% This default behavior can be modified by providing optional inputs.
%  
% =============================================================================
% =                                  Outputs                                  =
% =============================================================================
% Required Outputs:
%   x               Reconstructed signal.  For this documentation, we assume
%                   x has n total elements.  That is, it is of size compatable
%                   with the given A matrix/function call.  
% 
% Optional Outputs:
%   The optional outputs are in the following order:
%       optionalOutputs = [iter, objective, reconerror, cputime, solutionpath]
%
%   iter            The total number of iterations performed by the 
%                   algorithm.  Clearly this number will be between miniter
%                   and maxiter and will depend on the chosen stopping
%                   criteria. 
%
%   objective       The evolution of the objective function with the number
%                   of iterations.  The initial value of the objective
%                   function is stored in objective(1), and hence the 
%                   length of objective will be iter + 1.
%
%   reconerror      The evolution of the specified error metric with the
%                   number of iterations.  The reconstruction error can
%                   only be computed if the true underlying signal or image
%                   is provided using the 'TRUTH' option.  The error
%                   corresponding to the initial value is stored in
%                   reconerror(1), and hence the length of reconerror will
%                   be iter + 1.
%                   
%   cputime         Keeps track of the total elapsed time to reach each
%                   iteration.  This provides a better measure of the
%                   computational cost of the algorithm versus counting
%                   the number of iterations.  The clock starts at time
%                   cputime(1) = 0 and hence the length of cputime will
%                   also be iter + 1.
%
%   solutionpath    Provides a record of the intermediate iterates reached
%                   while computing the solution x.  Both the "noisy" step
%                   solutionpath.step and the "denoised" iterate
%                   solutionpath.iterate are saved.  The initialialization
%                   for the algorithm is stored in solutionpath(1).iterate.
%                   Since there is no corresponding initial step, 
%                   solutionpath(1).step contains all zeros.  Like all 
%                   the other output variables, the length of solutionpath
%                   will be iter + 1.
%


verbose = 0;
converged = 0;
iter = 1;  % Note: everything indexed by iter gives the value prior to the 
           % iterith iteration, e.g., objective(1) gives the objective
           % compted with xinit as the iterate
AT = [];
truth = [];
initialization = [];
warnings = 1;
recenter = 0;
mu = 0;
proj_2d =1;       % 1 do 2d projectios, 0 do 3d projections 
% Add a path to the denoising methods folder
spiraltapdir = which('SPIRALTAP');
[spiraltapdir dummy] = fileparts(spiraltapdir);
path([spiraltapdir,'/denoise'],path)

% ---- Noise Type ----
noisetype = 'Poisson';
% ---- For Poisson Noise ----
logepsilon = 1e-10;
sqrty = [];

% ---- Penalization Scheme ----
penalty = 'Canonical';
% l1 in an ONB
W = [];
WT = [];
subminiter = 1;
submaxiter = 50;
substopcriterion = 0;
subtolerance = 1e-5;
% Don't forget convergence criterion

% ---- For Choosing Alpha ----
alphamethod = 1;
monotone = 1;

% ---- Barz-Bor Scheme ---
alphainit = 1;
alphamin = 1e-30;
alphamax = 1e30;

% ---- For acceptance criterion ---
acceptalphamax = alphamax;
acceptmult = 2;
acceptdecrease = 0.1;
acceptpast = 10;

% ---- For termination ` ----
stopcriterion = 1;
miniter = 5;
maxiter = 100;
tolerance = 1e-6;

% ---- For Outputs ----
% By default, compute and save as little as possible:
saveobjective = 0;
computereconerror = 0; % by default assume f is not given
reconerrortype = 0; 
savecputime = 0;
savesolutionpath = 0;

% ---- Read in the input parameters ----
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for ii = 1:2:(length(varargin)-1)
        switch lower(varargin{ii})
            case 'verbose';             verbose             = varargin{ii+1};
            case 'at';                  AT                  = varargin{ii+1}; %
            case 'truth';               truth               = varargin{ii+1}; %
            case 'initialization';      initialization      = varargin{ii+1};
            case 'noisetype';           noisetype           = varargin{ii+1}; %
            case 'logepsilon';          logepsilon          = varargin{ii+1}; 
            case 'penalty';             penalty             = varargin{ii+1}; %
            case 'w';                   W                   = varargin{ii+1}; %
            case 'wt';                  WT                  = varargin{ii+1}; %
            case 'subminiter';          subminiter          = varargin{ii+1}; %
            case 'submaxiter';          submaxiter          = varargin{ii+1}; %
            case 'substopcriterion';    substopcriterion    = varargin{ii+1};
            case 'subtolerance';        subtolerance        = varargin{ii+1}; %
            case 'alphamethod';         alphamethod         = varargin{ii+1};
            case 'monotone';            monotone            = varargin{ii+1};
            case 'alphainit';           alphainit           = varargin{ii+1};
            case 'alphamin';            alphamin            = varargin{ii+1};
            case 'alphamax';            alphamax            = varargin{ii+1};
            case 'alphaaccept';         acceptalphamax      = varargin{ii+1};
            case 'eta';                 acceptmult          = varargin{ii+1};
            case 'sigma';               acceptdecrease      = varargin{ii+1};
            case 'acceptpast';          acceptpast          = varargin{ii+1};
            case 'stopcriterion';       stopcriterion       = varargin{ii+1};   
            case 'maxiter';             maxiter             = varargin{ii+1}; %
            case 'miniter';             miniter             = varargin{ii+1}; %
            case 'tolerance';           tolerance           = varargin{ii+1}; %
            case 'saveobjective';       saveobjective       = varargin{ii+1}; %
            case 'savereconerror';      savereconerror      = varargin{ii+1}; %
            case 'savecputime';         savecputime         = varargin{ii+1}; %
            case 'reconerrortype';      reconerrortype      = varargin{ii+1};
            case 'savesolutionpath';    savesolutionpath    = varargin{ii+1}; %
            case 'warnings';            warnings            = varargin{ii+1};
            case 'recenter';            recenter            = varargin{ii+1};
        otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''', varargin{ii}, '''']);
        end
    end
end

% ---- check the validity of the input parameters ----
% NOISETYPE:  For now only two options are available 'Poisson' and 'Gaussian'.
if sum( strcmpi(noisetype,{'poisson','gaussian'})) == 0
    error(['Invalid setting ''NOISETYPE'' = ''',num2str(noisetype),'''.  ',...
        'The parameter ''NOISETYPE'' may only be ''Gaussian'' or ''Poisson''.'])
end
% PENALTY:  The implemented penalty options are 'Canonical, 'ONB', 'RDP', 
% 'RDP-TI','TV'.
if sum( strcmpi(penalty,{'canonical','onb','rdp','rdp-ti','tv'})) == 0
    error(['Invalid setting ''PENALTY'' = ''',num2str(penalty),'''.  ',...
        'The parameter ''PENALTY'' may only be ''Canonical'', ''ONB'', ',...
        '''RDP'', ''RDP-TI'', or ''TV''.']);
end
% VERBOSE:  Needs to be a nonnegative integer.
if (round(verbose) ~= verbose) || (verbose < 0)
    error(['The parameter ''VERBOSE'' is required to be a nonnegative ',...
        'integer.  The setting ''VERBOSE'' = ',num2str(verbose),' is invalid.']);
end
% LOGEPSILON:  Needs to be nonnegative, usually small but that's relative.
if logepsilon < 0;
    error(['The parameter ''LOGEPSILON'' is required to be nonnegative.  ',...
        'The setting ''LOGEPSILON'' = ',num2str(tolerance),' is invalid.'])
end
% TOLERANCE:  Needs to be nonnegative, usually small but that's relative.
if tolerance < 0;
    error(['The parameter ''TOLERANCE'' is required to be nonnegative.  ',...
        'The setting ''TOLERANCE'' = ',num2str(tolerance),' is invalid.'])
end
% SUBTOLERANCE:  Needs to be nonnegative, usually small but that's relative.
if subtolerance < 0;
    error(['The parameter ''SUBTOLERANCE'' is required to be nonnegative.  ',...
        'The setting ''SUBTOLERANCE'' = ',num2str(subtolerance),' is invalid.'])
end
% MINITER and MAXITER:  Need to check that they are nonnegative integers and
% that miniter <= maxiter todo
if miniter > maxiter
    error(['The minimum number of iterations ''MINITER'' = ',...
        num2str(miniter),' exceeds the maximum number of iterations ',...
        '''MAXITER'' = ',num2str(maxiter),'.'])
end
if subminiter > submaxiter
    error(['The minimum number of subproblem iterations ''SUBMINITER'' = ',...
        num2str(subminiter),' exceeds the maximum number of subproblem ',...
        'iterations ''SUBMAXITER'' = ',num2str(submaxiter),'.'])
end
% AT:  If A is a matrix, AT is not required, but may optionally be provided.
% If A is a function call, AT is required.  In all cases, check that A and AT
% are of compatable size.  When A (and potentially AT) are given
% as matrices, we convert them to function calls for the remainder of the code
% Note: I think that it suffices to check whether or not the quantity
% dummy = y + A(AT(y)) is able to be computed, since it checks both the
% inner and outer dimensions of A and AT against that of the data y
if isa(A, 'function_handle') % A is a function call, so AT is required
    if isempty(AT) % AT simply not provided
        error(['Parameter ''AT'' not specified.  Please provide a method ',...
            'to compute A''*x matrix-vector products.'])
    else % AT was provided
        if isa(AT, 'function_handle') % A and AT are function calls
            try dummy = y + A(AT(y));
            catch exception; 
                error('Size incompatability between ''A'' and ''AT''.')
            end
        else % A is a function call, AT is a matrix        
            try dummy = y + A(AT*y);
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            AT = @(x) AT*x; % Define AT as a function call
        end
    end
else
    if isempty(AT) % A is a matrix, and AT not provided.
        AT = @(x) A'*x; % Just define function calls.
        A = @(x) A*x;
    else % A is a matrix, and AT provided, we need to check
        if isa(AT, 'function_handle') % A is a matrix, AT is a function call            
            try dummy = y + A*AT(y);
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            A = @(x) A*x; % Define A as a function call
        else % A and AT are matrices
            try dummy = y + A*AT*y;
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            AT = @(x) AT*x; % Define A and AT as function calls
            A = @(x) A*x;
        end
    end
end
% TRUTH:  Ensure that the size of truth, if given, is compatable with A and
% that it is nonnegative.  Note that this is irrespective of the noisetype
% since in the Gaussian case we still model the underlying signal as a
% nonnegative intensity.
if ~isempty(truth)
    try dummy = truth + AT(y);
    catch exception
        error(['The size of ''TRUTH'' is incompatable with the given ',...
            'sensing matrix ''A''.']);
    end
    if (min(truth(:)) < 0)
        error('The true signal ''TRUTH'' must be a nonnegative intensity.')
    end
end
% SAVEOBJECTIVE:  Just a binary indicator, check if not equal to 0 or 1.
if (numel(saveobjective) ~= 1)  || (sum( saveobjective == [0 1] ) ~= 1)
    error(['The option to save the objective evolution ',...
        'SAVEOBJECTIVE'' ',...
        'must be a binary scalar (either 0 or 1).'])
end     
% SAVERECONERROR:  Just a binary indicator, check if not equal to 0 or 1.
% If equal to 1, truth must be provided.
if (numel(savereconerror) ~= 1)  || (sum( savereconerror == [0 1] ) ~= 1)
    error(['The option to save the reconstruction error ',...
        'SAVERECONERROR'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
if savesolutionpath && isempty(truth)
    error(['The option to save the reconstruction error ',...
        '''SAVERECONERROR'' can only be used if the true signal ',...
        '''TRUTH'' is provided.'])
end
% SAVECPUTIME: Just a binary indicator, check if not equal to 0 or 1.
if (numel(savecputime) ~= 1)  || (sum( savecputime == [0 1] ) ~= 1)
    error(['The option to save the computation time ',...
        'SAVECPUTIME'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
% SAVESOLUTIONPATH: Just a binary indicator, check if not equal to 0 or 1.
if (numel(savesolutionpath) ~= 1)  || (sum( savesolutionpath == [0 1] ) ~= 1)
    error(['The option to save the solution path ',...
        'SAVESOLUTIONPATH'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
    

% Things to check and compute that depend on NOISETYPE:
switch lower(noisetype)
    case 'poisson'
        % Ensure that y is a vector of nonnegative counts
        if sum(round(y(:)) ~= y(:)) || (min(y(:)) < 0)
            error(['The data ''Y'' must contain nonnegative integer ',...
                'counts when ''NOISETYPE'' = ''Poisson''']);
        end
        % Maybe in future could check to ensure A and AT contain nonnegative
        %   elements, but perhaps too computationally wasteful 
        % Precompute useful quantities:
        sqrty = sqrt(y);
        % Ensure that recentering is not set
        if recenter
            todo
        end
    case 'gaussian'
        
end
% Things to check and compute that depend on PENALTY:
switch lower(penalty)
    case 'canonical'
        
    case 'onb' 
        % Already checked for valid subminiter, submaxiter, and subtolerance
        % Check for valid substopcriterion 
        % Need to check for the presense of W and WT
        if isempty(W)
            error(['Parameter ''W'' not specified.  Please provide a ',...
                'method to compute W*x matrix-vector products.'])
        end
        % Further checks to ensure we have both W and WT defined and that
        % the sizes are compatable by checking if y + A(WT(W(AT(y)))) can
        % be computed
        if isa(W, 'function_handle') % W is a function call, so WT is required
            if isempty(WT) % WT simply not provided
                error(['Parameter ''WT'' not specified.  Please provide a ',...
                    'method to compute W''*x matrix-vector products.'])
            else % WT was provided
        if isa(WT, 'function_handle') % W and WT are function calls
            try dummy = y + A(WT(W(AT(y))));
            catch exception; 
                error('Size incompatability between ''W'' and ''WT''.')
            end
        else % W is a function call, WT is a matrix        
            try dummy = y + A(WT*W(AT(y)));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            WT = @(x) WT*x; % Define WT as a function call
        end
    end
else
    if isempty(WT) % W is a matrix, and WT not provided.
        AT = @(x) W'*x; % Just define function calls.
        A = @(x) W*x;
    else % W is a matrix, and WT provided, we need to check
        if isa(WT, 'function_handle') % W is a matrix, WT is a function call            
            try dummy = y + A(WT(W*AT(y)));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            W = @(x) W*x; % Define W as a function call
        else % W and WT are matrices
            try dummy = y + A(WT(W*(AT(y))));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            WT = @(x) WT*x; % Define A and AT as function calls
            W = @(x) W*x;
        end
    end
end

	case 'rdp'
        %todo
        % Cannot enforce monotonicity (yet)
        if monotone
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Therefore monotonicity cannot be enforced.  ',...
                'Invalid option ''MONOTONIC'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        % Cannot compute objective function (yet)
        if saveobjective
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Invalid option ''SAVEOBJECTIVE'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
                
    case 'rdp-ti'
        % Cannot enforce monotonicity
        if monotone
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Therefore monotonicity cannot be enforced.  ',...
                'Invalid option ''MONOTONIC'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        % Cannot compute objective function 
        if saveobjective
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP-TI penalty.  ',...
                'Invalid option ''SAVEOBJECTIVE'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        
    case 'tv'
        % Cannot have a vectorized tau (yet)
        if (numel(tau) ~= 1)
            error(['A vector regularization parameter ''TAU'' cannot be ',...
                'used in conjuction with the TV penalty.']);
        end
end




% check that initialization is a scalar or a vector
% set initialization
if isempty(initialization);
    xinit = AT(y);
else
    xinit = initialization;
end

if recenter
    Aones = A(ones(size(xinit)));
    meanAones = mean(Aones(:));
    meany = mean(y(:));
    y = y - meany;
    mu = meany./meanAones;
    % Define new function calls for 'recentered' matrix
    A = @(x) A(x) - meanAones*sum(x(:))./length(xinit(:));
    AT = @(x) AT(x) - meanAones*sum(x(:))./length(xinit(:));
    % Adjust Initialization
    xinit = xinit - mu;
end
  

% ---- Check for validity of output parameters ----
% Check if there are too many or not enough
if (nargout == 0) && warnings
	disp('Warning:  You should reconsider not saving the output!');
	pause(1);
end
if (nargout < (2 + saveobjective + savereconerror ...
        + savecputime + savesolutionpath)) && warnings
    disp(['Warning:  Insufficient output parameters given to save ',...
        'the full output with the given options.']);
end
if nargout > (2 + saveobjective + savereconerror ...
        + savecputime + savesolutionpath)
        error('Too many output arguments specified for the given options.')
end

% --- Prepare for running the algorithm ----
% --- The below assumes that all parameters above are valid ---
% Initialize Main Algorithm 
x = xinit;
Ax = A(x);
alpha = alphainit;
Axprevious = Ax;
xprevious = x;
grad = computegrad(y,Ax,AT,noisetype,logepsilon);

% Prealocate arrays for storing results
% Initialize cputime and objective empty anyway (avoids errors in subfunctions):
cputime = [];
objective = [];

if savecputime
    cputime = zeros(maxiter+1,1);
end
if saveobjective
    objective = zeros(maxiter+1,1);
    objective(iter) = computeobjective(x,y,Ax,tau,noisetype,logepsilon,penalty,WT);
end
if savereconerror
    reconerror = zeros(maxiter+1,1);
    switch reconerrortype
        case 0 % RMS Error
            normtrue = sqrt( sum(truth(:).^2) );
            computereconerror = @(x) sqrt( sum( (x(:) + mu - truth(:) ).^2))./normtrue;
        case 1 % Relative absolute error
            normtrue = sum( abs(truth(:)) );
            computereconerror = @(x) sum( abs (x(:) + mu - truth(:)) )./normtrue;
    end
    reconerror(iter) = computereconerror(xinit);
end

if savesolutionpath
    % Note solutionpath(1).step will always be zeros since having an 
    % 'initial' step does not make sense
    solutionpath(1:maxiter+1) = struct('step',zeros(size(xinit)),...
        'iterate',zeros(size(xinit)));
    solutionpath(1).iterate = xinit;
end

if (verbose > 0)
    thetime = fix(clock);
    fprintf(['=========================================================\n',...
        '= Beginning SPIRAL Reconstruction    @ %2d:%2d %02d/%02d/%4d =\n',...
        '=   Noisetype: %-8s         Penalty: %-9s      =\n',...
        '=   Tau:       %-10.5e      Maxiter: %-5d          =\n',...
        '=========================================================\n'],...
        thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),...
        noisetype,penalty,tau,maxiter)      
end

tic; % Start clock for calculating computation time.
% =============================
% = Begin Main Algorithm Loop =
% =============================
while (iter <= miniter) || ((iter <= maxiter) && not(converged))
    disp("------------------- 3D Iteration --------------------")
    disp(iter)

    % ---- Compute the next iterate based on the method of computing alpha ----
    switch alphamethod
%         case 0 % Constant alpha throughout all iterations.
%             % If convergence criteria requires it, compute dx or dobjective
%             dx = xprevious;
%             step = xprevious - grad./alpha;
%             x = computesubsolution(step,tau,alpha,penalty,mu,...
%                 W,WT,subminiter,submaxiter,substopcriterion,...
%                 subtolerance);
%             dx = x - dx;
%             Ax = A(x);            
%             
        case 1 % Barzilai-Borwein choice of alpha
            if monotone 
                % do acceptance criterion.
                past = (max(iter-1-acceptpast,0):iter-1) + 1;
                maxpastobjective = max(objective(past));
                accept = 0;
                while (accept == 0)
                    
                    % --- Compute the step, and perform Gaussian 
                    %     denoising subproblem ----
                    dx = xprevious;
                    step = xprevious - grad./alpha;
                    x = computesubsolution_3d(step,tau,alpha,penalty,mu,iter,...
                        W,WT,subminiter,submaxiter,substopcriterion,...
                        subtolerance);
                    dx = x - dx;
                    Adx = Axprevious;
                    Ax = A(x);
                    Adx = Ax - Adx;
                    normsqdx = sum( dx(:).^2 );
                    
                    % --- Compute the resulting objective 
                    objective(iter + 1) = computeobjective(x,y,Ax,tau,...
                        noisetype,logepsilon,penalty,WT);
                        
                    if ( objective(iter+1) <= (maxpastobjective ...
                            - acceptdecrease*alpha/2*normsqdx) ) ...
                            || (alpha >= acceptalphamax);
                        accept = 1;
                    end
                    acceptalpha = alpha;  % Keep value for displaying
                    alpha = acceptmult*alpha;
                end
          
                    
            end
    end
    
    % ---- Calculate Output Quantities ----
    if savecputime
        cputime(iter+1) = toc;
    end
    if savereconerror
        reconerror(iter+1) = computereconerror(x);
    end
    if savesolutionpath
        solutionpath(iter+1).step = step;
        solutionpath(iter+1).iterate = x;
    end

    % Needed for next iteration and also termination criteria
    grad = computegrad(y,Ax,AT,noisetype,logepsilon);

    converged = checkconvergence(iter,miniter,stopcriterion,tolerance,...
                        dx, x, cputime(iter+1), objective);


    % Display progress
    if ~mod(iter,verbose)
        fprintf('Iter: %3d',iter);
        fprintf(', ||dx||%%: %11.4e', 100*norm(dx(:))/norm(x(:)));
        fprintf(', Alph: %11.4e',alpha);
        if monotone
            fprintf(', Alph Acc: %11.4e',acceptalpha)
        end
        if savecputime
            fprintf(', Time: %3d',cputime(iter+1));
        end
        if saveobjective
            fprintf(', Obj: %11.4e',objective(iter+1));
            fprintf(', dObj%%: %11.4e',...
                100*abs(objective(iter+1) - objective(iter))./...
                abs(objective(iter)))
        end      
        if savereconerror
            fprintf(', Err: %11.4e',reconerror(iter+1))
        end
        fprintf('\n')
    end 
    
    
    % --- Prepare for next iteration ---- 
    % Update alpha
    switch alphamethod
        case 0 % do nothing, constant alpha
        case 1 % bb method
            %Adx is overwritten at top of iteration, so this is an ok reuse
            % Adx is overwritten at top of iteration, so this is an ok reuse
            switch lower(noisetype)
                case 'poisson'
                    Adx = Adx.*sqrty./(Ax + logepsilon); 
                case 'gaussian'
                    % No need to scale Adx
            end
            gamma = sum(Adx(:).^2);
            if gamma == 0
                alpha = alphamin;
            else
                alpha = gamma./normsqdx;
                alpha = min(alphamax, max(alpha, alphamin));
            end
    end
    
    % --- Store current values as previous values for next iteration ---
    xprevious = x;
    Axprevious = Ax; 
    iter = iter + 1;
end
% ===========================
% = End Main Algorithm Loop =
% ===========================
% Add on mean if recentered (if not mu == 0);
x = x + mu;

% Determine what needs to be in the variable output and
% crop the output if the maximum number of iterations were not used.
% Note, need to subtract 1 since iter is incremented at the end of the loop
iter = iter - 1;
varargout = {iter};

if saveobjective
    varargout = [varargout {objective(1:iter+1);}];
end
if savereconerror
    varargout = [varargout {reconerror(1:iter+1)}];
end
if savecputime
    varargout = [varargout {cputime(1:iter+1)}];
end
if savesolutionpath
    varargout = [varargout {solutionpath(1:iter+1)}];
end

if (verbose > 0)
    thetime = fix(clock);
    fprintf(['=========================================================\n',...
        '= Completed SPIRAL Reconstruction    @ %2d:%2d %02d/%02d/%4d =\n',...
        '=   Noisetype: %-8s         Penalty: %-9s      =\n',...
        '=   Tau:       %-10.5e      Iter:    %-5d          =\n'],...
        thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),...
        noisetype,penalty,tau,iter)      
    fprintf('=========================================================\n');
end

end

function [x, varargout] = SPIRALTAP_1GP1P1C_Hap_2d(y, A, tau, varargin)
% ==== Set default/initial parameter values ====
% ---- All Methods ----
% =============================================================================
% =        SPIRAL:  Sparse Poisson Intensity Reconstruction Algorithms        =
% =                                Version 1.0                                =
% =============================================================================
% =    Copyright 2009, 2010                                                   =
% =    Zachary T. Harmany*, Roummel F. Marcia**, Rebecca M. Willett*          =
% =        *  Department of Electrical and Computer Engineering               =
% =           Duke University                                                 =
% =           Durham, NC 27708, USA                                           =
% =       **  School of Natural Sciences                                      =
% =           University of California, Merced                                =
% =           Merced, CA 95343, USA                                           =
% =                                                                           =
% =    Corresponding author: Zachary T. Harmany (zth@duke.edu)                =
% ============================================================================= 
%
% =============================================================================
% =                               Documentation                               =
% =============================================================================
% Syntax:
%   [x, optionalOutputs] = SPIRALTAP(y, A, tau, optionalInputs)
% 
%   More details and supporting publications are 
%   available on the SPIRAL Toolbox homepage
%   http://www.ee.duke.edu/~zth/spiral/
% 
% =============================================================================
% =                                  Inputs                                   =
% =============================================================================
% Required Inputs:
%	y               Degraded observations.  For this documenation, we say
%                   that y has m total elements.
%
%   A               Sensing / observation matrix.  For this documentation, 
%                   we say that A is an m x n matrix.  A could also be a
%                   function call A() that computes matrix-vector products
%                   such that A(x) = A*x where x has n total elements.  In 
%                   this case one must also specify a function call AT() 
%                   using the 'AT' option that computes matrix-vector 
%                   products with the adjoint of A such that AT(x) = A'*x.  
%
%   tau             Regularization parameter that trades off the data fit
%                   (negative log-likelihood) with the regularization.
%                   The regularization parameter can either be a
%                   nonnegative real scalar or (for all methods except the
%                   total variation penalty) have n nonnegative real 
%                   elements which allows for nonuniform penalization 
%                   schemes. 
%           
% Optional Inputs:
% If one were to only input y, A, and tau into the algorithm, there are 
% many necessary assumptions as to what to do with the inputs.  By default
% SPIRAL assumes that:
%   - y contains Poisson realizations of A*f, where f is the true underlying
%     signal (to be estimated),
%   - the penalty is the l_1 norm of x (i.e., we promote sparsity in the 
%     canonical basis.
% This default behavior can be modified by providing optional inputs.
%  
% =============================================================================
% =                                  Outputs                                  =
% =============================================================================
% Required Outputs:
%   x               Reconstructed signal.  For this documentation, we assume
%                   x has n total elements.  That is, it is of size compatable
%                   with the given A matrix/function call.  
% 
% Optional Outputs:
%   The optional outputs are in the following order:
%       optionalOutputs = [iter, objective, reconerror, cputime, solutionpath]
%
%   iter            The total number of iterations performed by the 
%                   algorithm.  Clearly this number will be between miniter
%                   and maxiter and will depend on the chosen stopping
%                   criteria. 
%
%   objective       The evolution of the objective function with the number
%                   of iterations.  The initial value of the objective
%                   function is stored in objective(1), and hence the 
%                   length of objective will be iter + 1.
%
%   reconerror      The evolution of the specified error metric with the
%                   number of iterations.  The reconstruction error can
%                   only be computed if the true underlying signal or image
%                   is provided using the 'TRUTH' option.  The error
%                   corresponding to the initial value is stored in
%                   reconerror(1), and hence the length of reconerror will
%                   be iter + 1.
%                   
%   cputime         Keeps track of the total elapsed time to reach each
%                   iteration.  This provides a better measure of the
%                   computational cost of the algorithm versus counting
%                   the number of iterations.  The clock starts at time
%                   cputime(1) = 0 and hence the length of cputime will
%                   also be iter + 1.
%
%   solutionpath    Provides a record of the intermediate iterates reached
%                   while computing the solution x.  Both the "noisy" step
%                   solutionpath.step and the "denoised" iterate
%                   solutionpath.iterate are saved.  The initialialization
%                   for the algorithm is stored in solutionpath(1).iterate.
%                   Since there is no corresponding initial step, 
%                   solutionpath(1).step contains all zeros.  Like all 
%                   the other output variables, the length of solutionpath
%                   will be iter + 1.
%


verbose = 0;
converged = 0;
iter = 1;  % Note: everything indexed by iter gives the value prior to the 
           % iterith iteration, e.g., objective(1) gives the objective
           % compted with xinit as the iterate
AT = [];
truth = [];
initialization = [];
warnings = 1;
recenter = 0;
mu = 0;
proj_2d =1;       % 1 do 2d projectios, 0 do 3d projections 
% Add a path to the denoising methods folder
spiraltapdir = which('SPIRALTAP');
[spiraltapdir dummy] = fileparts(spiraltapdir);
path([spiraltapdir,'/denoise'],path)

% ---- Noise Type ----
noisetype = 'Poisson';
% ---- For Poisson Noise ----
logepsilon = 1e-10;
sqrty = [];

% ---- Penalization Scheme ----
penalty = 'Canonical';
% l1 in an ONB
W = [];
WT = [];
subminiter = 1;
submaxiter = 50;
substopcriterion = 0;
subtolerance = 1e-5;
% Don't forget convergence criterion

% ---- For Choosing Alpha ----
alphamethod = 1;
monotone = 1;

% ---- Barz-Bor Scheme ---
alphainit = 1;
alphamin = 1e-30;
alphamax = 1e30;

% ---- For acceptance criterion ---
acceptalphamax = alphamax;
acceptmult = 2;
acceptdecrease = 0.1;
acceptpast = 10;

% ---- For termination ` ----
stopcriterion = 1;
miniter = 5;
maxiter = 100;
tolerance = 1e-6;

% ---- For Outputs ----
% By default, compute and save as little as possible:
saveobjective = 0;
computereconerror = 0; % by default assume f is not given
reconerrortype = 0; 
savecputime = 0;
savesolutionpath = 0;

% ---- Read in the input parameters ----
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for ii = 1:2:(length(varargin)-1)
        switch lower(varargin{ii})
            case 'verbose';             verbose             = varargin{ii+1};
            case 'at';                  AT                  = varargin{ii+1}; %
            case 'truth';               truth               = varargin{ii+1}; %
            case 'initialization';      initialization      = varargin{ii+1};
            case 'noisetype';           noisetype           = varargin{ii+1}; %
            case 'logepsilon';          logepsilon          = varargin{ii+1}; 
            case 'penalty';             penalty             = varargin{ii+1}; %
            case 'w';                   W                   = varargin{ii+1}; %
            case 'wt';                  WT                  = varargin{ii+1}; %
            case 'subminiter';          subminiter          = varargin{ii+1}; %
            case 'submaxiter';          submaxiter          = varargin{ii+1}; %
            case 'substopcriterion';    substopcriterion    = varargin{ii+1};
            case 'subtolerance';        subtolerance        = varargin{ii+1}; %
            case 'alphamethod';         alphamethod         = varargin{ii+1};
            case 'monotone';            monotone            = varargin{ii+1};
            case 'alphainit';           alphainit           = varargin{ii+1};
            case 'alphamin';            alphamin            = varargin{ii+1};
            case 'alphamax';            alphamax            = varargin{ii+1};
            case 'alphaaccept';         acceptalphamax      = varargin{ii+1};
            case 'eta';                 acceptmult          = varargin{ii+1};
            case 'sigma';               acceptdecrease      = varargin{ii+1};
            case 'acceptpast';          acceptpast          = varargin{ii+1};
            case 'stopcriterion';       stopcriterion       = varargin{ii+1};   
            case 'maxiter';             maxiter             = varargin{ii+1}; %
            case 'miniter';             miniter             = varargin{ii+1}; %
            case 'tolerance';           tolerance           = varargin{ii+1}; %
            case 'saveobjective';       saveobjective       = varargin{ii+1}; %
            case 'savereconerror';      savereconerror      = varargin{ii+1}; %
            case 'savecputime';         savecputime         = varargin{ii+1}; %
            case 'reconerrortype';      reconerrortype      = varargin{ii+1};
            case 'savesolutionpath';    savesolutionpath    = varargin{ii+1}; %
            case 'warnings';            warnings            = varargin{ii+1};
            case 'recenter';            recenter            = varargin{ii+1};
        otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''', varargin{ii}, '''']);
        end
    end
end

% ---- check the validity of the input parameters ----
% NOISETYPE:  For now only two options are available 'Poisson' and 'Gaussian'.
if sum( strcmpi(noisetype,{'poisson','gaussian'})) == 0
    error(['Invalid setting ''NOISETYPE'' = ''',num2str(noisetype),'''.  ',...
        'The parameter ''NOISETYPE'' may only be ''Gaussian'' or ''Poisson''.'])
end
% PENALTY:  The implemented penalty options are 'Canonical, 'ONB', 'RDP', 
% 'RDP-TI','TV'.
if sum( strcmpi(penalty,{'canonical','onb','rdp','rdp-ti','tv'})) == 0
    error(['Invalid setting ''PENALTY'' = ''',num2str(penalty),'''.  ',...
        'The parameter ''PENALTY'' may only be ''Canonical'', ''ONB'', ',...
        '''RDP'', ''RDP-TI'', or ''TV''.']);
end
% VERBOSE:  Needs to be a nonnegative integer.
if (round(verbose) ~= verbose) || (verbose < 0)
    error(['The parameter ''VERBOSE'' is required to be a nonnegative ',...
        'integer.  The setting ''VERBOSE'' = ',num2str(verbose),' is invalid.']);
end
% LOGEPSILON:  Needs to be nonnegative, usually small but that's relative.
if logepsilon < 0;
    error(['The parameter ''LOGEPSILON'' is required to be nonnegative.  ',...
        'The setting ''LOGEPSILON'' = ',num2str(tolerance),' is invalid.'])
end
% TOLERANCE:  Needs to be nonnegative, usually small but that's relative.
if tolerance < 0;
    error(['The parameter ''TOLERANCE'' is required to be nonnegative.  ',...
        'The setting ''TOLERANCE'' = ',num2str(tolerance),' is invalid.'])
end
% SUBTOLERANCE:  Needs to be nonnegative, usually small but that's relative.
if subtolerance < 0;
    error(['The parameter ''SUBTOLERANCE'' is required to be nonnegative.  ',...
        'The setting ''SUBTOLERANCE'' = ',num2str(subtolerance),' is invalid.'])
end
% MINITER and MAXITER:  Need to check that they are nonnegative integers and
% that miniter <= maxiter todo
if miniter > maxiter
    error(['The minimum number of iterations ''MINITER'' = ',...
        num2str(miniter),' exceeds the maximum number of iterations ',...
        '''MAXITER'' = ',num2str(maxiter),'.'])
end
if subminiter > submaxiter
    error(['The minimum number of subproblem iterations ''SUBMINITER'' = ',...
        num2str(subminiter),' exceeds the maximum number of subproblem ',...
        'iterations ''SUBMAXITER'' = ',num2str(submaxiter),'.'])
end
% AT:  If A is a matrix, AT is not required, but may optionally be provided.
% If A is a function call, AT is required.  In all cases, check that A and AT
% are of compatable size.  When A (and potentially AT) are given
% as matrices, we convert them to function calls for the remainder of the code
% Note: I think that it suffices to check whether or not the quantity
% dummy = y + A(AT(y)) is able to be computed, since it checks both the
% inner and outer dimensions of A and AT against that of the data y
if isa(A, 'function_handle') % A is a function call, so AT is required
    if isempty(AT) % AT simply not provided
        error(['Parameter ''AT'' not specified.  Please provide a method ',...
            'to compute A''*x matrix-vector products.'])
    else % AT was provided
        if isa(AT, 'function_handle') % A and AT are function calls
            try dummy = y + A(AT(y));
            catch exception; 
                error('Size incompatability between ''A'' and ''AT''.')
            end
        else % A is a function call, AT is a matrix        
            try dummy = y + A(AT*y);
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            AT = @(x) AT*x; % Define AT as a function call
        end
    end
else
    if isempty(AT) % A is a matrix, and AT not provided.
        AT = @(x) A'*x; % Just define function calls.
        A = @(x) A*x;
    else % A is a matrix, and AT provided, we need to check
        if isa(AT, 'function_handle') % A is a matrix, AT is a function call            
            try dummy = y + A*AT(y);
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            A = @(x) A*x; % Define A as a function call
        else % A and AT are matrices
            try dummy = y + A*AT*y;
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            AT = @(x) AT*x; % Define A and AT as function calls
            A = @(x) A*x;
        end
    end
end
% TRUTH:  Ensure that the size of truth, if given, is compatable with A and
% that it is nonnegative.  Note that this is irrespective of the noisetype
% since in the Gaussian case we still model the underlying signal as a
% nonnegative intensity.
if ~isempty(truth)
    try dummy = truth + AT(y);
    catch exception
        error(['The size of ''TRUTH'' is incompatable with the given ',...
            'sensing matrix ''A''.']);
    end
    if (min(truth(:)) < 0)
        error('The true signal ''TRUTH'' must be a nonnegative intensity.')
    end
end
% SAVEOBJECTIVE:  Just a binary indicator, check if not equal to 0 or 1.
if (numel(saveobjective) ~= 1)  || (sum( saveobjective == [0 1] ) ~= 1)
    error(['The option to save the objective evolution ',...
        'SAVEOBJECTIVE'' ',...
        'must be a binary scalar (either 0 or 1).'])
end     
% SAVERECONERROR:  Just a binary indicator, check if not equal to 0 or 1.
% If equal to 1, truth must be provided.
if (numel(savereconerror) ~= 1)  || (sum( savereconerror == [0 1] ) ~= 1)
    error(['The option to save the reconstruction error ',...
        'SAVERECONERROR'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
if savesolutionpath && isempty(truth)
    error(['The option to save the reconstruction error ',...
        '''SAVERECONERROR'' can only be used if the true signal ',...
        '''TRUTH'' is provided.'])
end
% SAVECPUTIME: Just a binary indicator, check if not equal to 0 or 1.
if (numel(savecputime) ~= 1)  || (sum( savecputime == [0 1] ) ~= 1)
    error(['The option to save the computation time ',...
        'SAVECPUTIME'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
% SAVESOLUTIONPATH: Just a binary indicator, check if not equal to 0 or 1.
if (numel(savesolutionpath) ~= 1)  || (sum( savesolutionpath == [0 1] ) ~= 1)
    error(['The option to save the solution path ',...
        'SAVESOLUTIONPATH'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
    

% Things to check and compute that depend on NOISETYPE:
switch lower(noisetype)
    case 'poisson'
        % Ensure that y is a vector of nonnegative counts
        if sum(round(y(:)) ~= y(:)) || (min(y(:)) < 0)
            error(['The data ''Y'' must contain nonnegative integer ',...
                'counts when ''NOISETYPE'' = ''Poisson''']);
        end
        % Maybe in future could check to ensure A and AT contain nonnegative
        %   elements, but perhaps too computationally wasteful 
        % Precompute useful quantities:
        sqrty = sqrt(y);
        % Ensure that recentering is not set
        if recenter
            todo
        end
    case 'gaussian'
        
end
% Things to check and compute that depend on PENALTY:
switch lower(penalty)
    case 'canonical'
        
    case 'onb' 
        % Already checked for valid subminiter, submaxiter, and subtolerance
        % Check for valid substopcriterion 
        % Need to check for the presense of W and WT
        if isempty(W)
            error(['Parameter ''W'' not specified.  Please provide a ',...
                'method to compute W*x matrix-vector products.'])
        end
        % Further checks to ensure we have both W and WT defined and that
        % the sizes are compatable by checking if y + A(WT(W(AT(y)))) can
        % be computed
        if isa(W, 'function_handle') % W is a function call, so WT is required
            if isempty(WT) % WT simply not provided
                error(['Parameter ''WT'' not specified.  Please provide a ',...
                    'method to compute W''*x matrix-vector products.'])
            else % WT was provided
        if isa(WT, 'function_handle') % W and WT are function calls
            try dummy = y + A(WT(W(AT(y))));
            catch exception; 
                error('Size incompatability between ''W'' and ''WT''.')
            end
        else % W is a function call, WT is a matrix        
            try dummy = y + A(WT*W(AT(y)));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            WT = @(x) WT*x; % Define WT as a function call
        end
    end
else
    if isempty(WT) % W is a matrix, and WT not provided.
        AT = @(x) W'*x; % Just define function calls.
        A = @(x) W*x;
    else % W is a matrix, and WT provided, we need to check
        if isa(WT, 'function_handle') % W is a matrix, WT is a function call            
            try dummy = y + A(WT(W*AT(y)));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            W = @(x) W*x; % Define W as a function call
        else % W and WT are matrices
            try dummy = y + A(WT(W*(AT(y))));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            WT = @(x) WT*x; % Define A and AT as function calls
            W = @(x) W*x;
        end
    end
end

	case 'rdp'
        %todo
        % Cannot enforce monotonicity (yet)
        if monotone
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Therefore monotonicity cannot be enforced.  ',...
                'Invalid option ''MONOTONIC'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        % Cannot compute objective function (yet)
        if saveobjective
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Invalid option ''SAVEOBJECTIVE'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
                
    case 'rdp-ti'
        % Cannot enforce monotonicity
        if monotone
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Therefore monotonicity cannot be enforced.  ',...
                'Invalid option ''MONOTONIC'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        % Cannot compute objective function 
        if saveobjective
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP-TI penalty.  ',...
                'Invalid option ''SAVEOBJECTIVE'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        
    case 'tv'
        % Cannot have a vectorized tau (yet)
        if (numel(tau) ~= 1)
            error(['A vector regularization parameter ''TAU'' cannot be ',...
                'used in conjuction with the TV penalty.']);
        end
end




% check that initialization is a scalar or a vector
% set initialization
if isempty(initialization);
    xinit = AT(y);
else
    xinit = initialization;
end

if recenter
    Aones = A(ones(size(xinit)));
    meanAones = mean(Aones(:));
    meany = mean(y(:));
    y = y - meany;
    mu = meany./meanAones;
    % Define new function calls for 'recentered' matrix
    A = @(x) A(x) - meanAones*sum(x(:))./length(xinit(:));
    AT = @(x) AT(x) - meanAones*sum(x(:))./length(xinit(:));
    % Adjust Initialization
    xinit = xinit - mu;
end
  

% ---- Check for validity of output parameters ----
% Check if there are too many or not enough
if (nargout == 0) && warnings
	disp('Warning:  You should reconsider not saving the output!');
	pause(1);
end
if (nargout < (2 + saveobjective + savereconerror ...
        + savecputime + savesolutionpath)) && warnings
    disp(['Warning:  Insufficient output parameters given to save ',...
        'the full output with the given options.']);
end
if nargout > (2 + saveobjective + savereconerror ...
        + savecputime + savesolutionpath)
        error('Too many output arguments specified for the given options.')
end

% --- Prepare for running the algorithm ----
% --- The below assumes that all parameters above are valid ---
% Initialize Main Algorithm 
x = xinit;
Ax = A(x);
alpha = alphainit;
Axprevious = Ax;
xprevious = x;
grad = computegrad(y,Ax,AT,noisetype,logepsilon);

% Prealocate arrays for storing results
% Initialize cputime and objective empty anyway (avoids errors in subfunctions):
cputime = [];
objective = [];

if savecputime
    cputime = zeros(maxiter+1,1);
end
if saveobjective
    objective = zeros(maxiter+1,1);
    objective(iter) = computeobjective(x,y,Ax,tau,noisetype,logepsilon,penalty,WT);
end
if savereconerror
    reconerror = zeros(maxiter+1,1);
    switch reconerrortype
        case 0 % RMS Error
            normtrue = sqrt( sum(truth(:).^2) );
            computereconerror = @(x) sqrt( sum( (x(:) + mu - truth(:) ).^2))./normtrue;
        case 1 % Relative absolute error
            normtrue = sum( abs(truth(:)) );
            computereconerror = @(x) sum( abs (x(:) + mu - truth(:)) )./normtrue;
    end
    reconerror(iter) = computereconerror(xinit);
end

if savesolutionpath
    % Note solutionpath(1).step will always be zeros since having an 
    % 'initial' step does not make sense
    solutionpath(1:maxiter+1) = struct('step',zeros(size(xinit)),...
        'iterate',zeros(size(xinit)));
    solutionpath(1).iterate = xinit;
end

if (verbose > 0)
    thetime = fix(clock);
    fprintf(['=========================================================\n',...
        '= Beginning SPIRAL Reconstruction    @ %2d:%2d %02d/%02d/%4d =\n',...
        '=   Noisetype: %-8s         Penalty: %-9s      =\n',...
        '=   Tau:       %-10.5e      Maxiter: %-5d          =\n',...
        '=========================================================\n'],...
        thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),...
        noisetype,penalty,tau,maxiter)      
end

tic; % Start clock for calculating computation time.
% =============================
% = Begin Main Algorithm Loop =
% =============================
while (iter <= miniter) || ((iter <= maxiter) && not(converged))

    disp("================")
    disp(iter)
    % ---- Compute the next iterate based on the method of computing alpha ----
    switch alphamethod
%         case 0 % Constant alpha throughout all iterations.
%             % If convergence criteria requires it, compute dx or dobjective
%             dx = xprevious;
%             step = xprevious - grad./alpha;
%             x = computesubsolution(step,tau,alpha,penalty,mu,...
%                 W,WT,subminiter,submaxiter,substopcriterion,...
%                 subtolerance);
%             dx = x - dx;
%             Ax = A(x);            
%             
        case 1 % Barzilai-Borwein choice of alpha
            if monotone 
                % do acceptance criterion.
                past = (max(iter-1-acceptpast,0):iter-1) + 1;
                maxpastobjective = max(objective(past));
                accept = 0;
                while (accept == 0)
                    
                    % --- Compute the step, and perform Gaussian 
                    %     denoising subproblem ----
                    dx = xprevious;
                    step = xprevious - grad./alpha;
                    x = computesubsolution_2d(step,tau,alpha,penalty,mu,iter,...
                        W,WT,subminiter,submaxiter,substopcriterion,...
                        subtolerance);
                    dx = x - dx;
                    Adx = Axprevious;
                    Ax = A(x);
                    Adx = Ax - Adx;
                    normsqdx = sum( dx(:).^2 );
                    
                    % --- Compute the resulting objective 
                    objective(iter + 1) = computeobjective(x,y,Ax,tau,...
                        noisetype,logepsilon,penalty,WT);
                        
                    if ( objective(iter+1) <= (maxpastobjective ...
                            - acceptdecrease*alpha/2*normsqdx) ) ...
                            || (alpha >= acceptalphamax);
                        accept = 1;
                    end
                    acceptalpha = alpha;  % Keep value for displaying
                    alpha = acceptmult*alpha;
                end
          
                    
            end
    end
    
    % ---- Calculate Output Quantities ----
    if savecputime
        cputime(iter+1) = toc;
    end
    if savereconerror
        reconerror(iter+1) = computereconerror(x);
    end
    if savesolutionpath
        solutionpath(iter+1).step = step;
        solutionpath(iter+1).iterate = x;
    end

    % Needed for next iteration and also termination criteria
    grad = computegrad(y,Ax,AT,noisetype,logepsilon);

    converged = checkconvergence(iter,miniter,stopcriterion,tolerance,...
                        dx, x, cputime(iter+1), objective);


    % Display progress
    if ~mod(iter,verbose)
        fprintf('Iter: %3d',iter);
        fprintf(', ||dx||%%: %11.4e', 100*norm(dx(:))/norm(x(:)));
        fprintf(', Alph: %11.4e',alpha);
        if monotone
            fprintf(', Alph Acc: %11.4e',acceptalpha)
        end
        if savecputime
            fprintf(', Time: %3d',cputime(iter+1));
        end
        if saveobjective
            fprintf(', Obj: %11.4e',objective(iter+1));
            fprintf(', dObj%%: %11.4e',...
                100*abs(objective(iter+1) - objective(iter))./...
                abs(objective(iter)))
        end      
        if savereconerror
            fprintf(', Err: %11.4e',reconerror(iter+1))
        end
        fprintf('\n')
    end 
    
    
    % --- Prepare for next iteration ---- 
    % Update alpha
    switch alphamethod
        case 0 % do nothing, constant alpha
        case 1 % bb method
            %Adx is overwritten at top of iteration, so this is an ok reuse
            % Adx is overwritten at top of iteration, so this is an ok reuse
            switch lower(noisetype)
                case 'poisson'
                    Adx = Adx.*sqrty./(Ax + logepsilon); 
                case 'gaussian'
                    % No need to scale Adx
            end
            gamma = sum(Adx(:).^2);
            if gamma == 0
                alpha = alphamin;
            else
                alpha = gamma./normsqdx;
                alpha = min(alphamax, max(alpha, alphamin));
            end
    end
    
    % --- Store current values as previous values for next iteration ---
    xprevious = x;
    Axprevious = Ax; 
    iter = iter + 1;
end
% ===========================
% = End Main Algorithm Loop =
% ===========================
% Add on mean if recentered (if not mu == 0);
x = x + mu;

% Determine what needs to be in the variable output and
% crop the output if the maximum number of iterations were not used.
% Note, need to subtract 1 since iter is incremented at the end of the loop
iter = iter - 1;
varargout = {iter};

if saveobjective
    varargout = [varargout {objective(1:iter+1);}];
end
if savereconerror
    varargout = [varargout {reconerror(1:iter+1)}];
end
if savecputime
    varargout = [varargout {cputime(1:iter+1)}];
end
if savesolutionpath
    varargout = [varargout {solutionpath(1:iter+1)}];
end

if (verbose > 0)
    thetime = fix(clock);
    fprintf(['=========================================================\n',...
        '= Completed SPIRAL Reconstruction    @ %2d:%2d %02d/%02d/%4d =\n',...
        '=   Noisetype: %-8s         Penalty: %-9s      =\n',...
        '=   Tau:       %-10.5e      Iter:    %-5d          =\n'],...
        thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),...
        noisetype,penalty,tau,iter)      
    fprintf('=========================================================\n');
end

end

% =============================================================================
% =============================================================================
% =============================================================================
% =                            Helper Subfunctions                            =
% =============================================================================
% =============================================================================
% =============================================================================

% ==========================
% = Gradient Computation: =
% =========================
function grad = computegrad(y,Ax,AT,noisetype,logepsilon)
    % Perhaps change to varargin 
    switch lower(noisetype)
        case 'poisson'
            grad = AT(1 - (y./(Ax + logepsilon)));
        case 'gaussian'
            grad = AT(Ax - y);
    end
end

% ==========================
% = Objective Computation: =
% ==========================
function objective = computeobjective(x,y,Ax,tau,noisetype,logepsilon,...
    penalty,varargin)
% Perhaps change to varargin 
% 1) Compute log-likelihood:
switch lower(noisetype)
    case 'poisson'
        precompute = y.*log(Ax + logepsilon);
        objective = sum(Ax(:)) - sum(precompute(:));
    case 'gaussian'
        objective = sum( (y(:) - Ax(:)).^2)./2;
end
% 2) Compute Penalty:
switch lower(penalty)
    case 'canonical'
        n = length(x)/3;
        objective = objective + sum(abs(tau(1).*x(1:n))); % use tau(1) for child inherited variants
        objective = objective + sum(abs(tau(2).*x(n+1:2*n))); % use tau(2) for parent variants
        objective = objective + sum(abs(tau(3).*x(2*n+1:3*n))); % use tau(3) for grand parent variants
%        objective = objective + sum(abs(tau(:).*x(:)));
	case 'onb' 
    	WT = varargin{1};
        WTx = WT(x);
        objective = objective + sum(abs(tau(:).*WTx(:)));
	case 'rdp'
        todo
    case 'rdp-ti'
        todo
    case 'tv'
        objective = objective + tau.*tlv(x,'l1');
end
end

% =====================================
% = Denoising Subproblem Computation: =
% =====================================
function subsolution = computesubsolution_3d(step,tau,alpha,penalty,mu,iter,varargin)
    switch lower(penalty)
        case 'canonical'
            
            %%%subsolution = max(step - tau./alpha + mu, 0.0);
            %%%subsolution = max(step - tau./alpha + mu, 0.0);
            %subsolution = step - tau./alpha;
            n = length(step)/3;
            subsolution_c = step(1:n) - tau(1)./alpha; % child
            % variants inherited from parent
            subsolution_p = step(n+1:2*n) - tau(2)./alpha; % novel
            % child variants
            subsolution_gp = step(2*n+1:3*n) - tau(3)./alpha; % parent
            % variants
            %subsolution = Novel_const(subsolution);
            subsolution=[subsolution_c;subsolution_p;subsolution_gp];
            
            %if proj_2d == 1   %2D projections 
             %[subsolution] = Projections_2D_Alternating(subsolution,10);      
            %else % do 3D projections 
             %[subsolution] = Projct_3D(subsolution);                   
            %end 
             
        case 'onb'
            % if onb is selected, varargin must be such that
            W                   = varargin{1};
            WT                  = varargin{2};
            subminiter          = varargin{3};
            submaxiter          = varargin{4};
            substopcriterion    = varargin{5};
            subtolerance        = varargin{6};
                                   
            subsolution = constrainedl2l1denoise(step,W,WT,tau./alpha,mu,...
                subminiter,submaxiter,substopcriterion,subtolerance);
        case 'rdp'
            subsolution = haarTVApprox2DNN_recentered(step,tau./alpha,-mu);
        case 'rdp-ti'
            subsolution = haarTIApprox2DNN_recentered(step,tau./alpha,-mu);
        case 'tv'
            subtolerance        = varargin{6};
            submaxiter          = varargin{4};
            % From Becca's Code:
            pars.print = 0;
            pars.tv = 'l1';
            pars.MAXITER = submaxiter;
            pars.epsilon = subtolerance; % Becca used 1e-5;
            if tau>0
                subsolution = denoise_bound(step,tau./alpha,-mu,Inf,pars);
            else
                subsolution = step.*(step>0);
            end
    end           
end
function subsolution = computesubsolution_2d(step,tau,alpha,penalty,mu,iter,varargin)
    switch lower(penalty)
        case 'canonical'
            
            %%%subsolution = max(step - tau./alpha + mu, 0.0);
            %%%subsolution = max(step - tau./alpha + mu, 0.0);
            %subsolution = step - tau./alpha;
            n = length(step)/3;
            subsolution_c = step(1:n) - tau(1)./alpha; % child
            % variants inherited from parent
            subsolution_p = step(n+1:2*n) - tau(2)./alpha; % novel
            % child variants
            subsolution_gp = step(2*n+1:3*n) - tau(3)./alpha; % parent
            % variants
            %subsolution = Novel_const(subsolution);
            subsolution=[subsolution_c;subsolution_p;subsolution_gp];
            
            %if proj_2d == 1   %2D projections 
             %[subsolution] = Alt2D_Proj(subsolution,10);      
            %else % do 3D projections 
             %[subsolution] = Projections_3D(subsolution);                   
            %end 
             
        case 'onb'
            % if onb is selected, varargin must be such that
            W                   = varargin{1};
            WT                  = varargin{2};
            subminiter          = varargin{3};
            submaxiter          = varargin{4};
            substopcriterion    = varargin{5};
            subtolerance        = varargin{6};
                                   
            subsolution = constrainedl2l1denoise(step,W,WT,tau./alpha,mu,...
                subminiter,submaxiter,substopcriterion,subtolerance);
        case 'rdp'
            subsolution = haarTVApprox2DNN_recentered(step,tau./alpha,-mu);
        case 'rdp-ti'
            subsolution = haarTIApprox2DNN_recentered(step,tau./alpha,-mu);
        case 'tv'
            subtolerance        = varargin{6};
            submaxiter          = varargin{4};
            % From Becca's Code:
            pars.print = 0;
            pars.tv = 'l1';
            pars.MAXITER = submaxiter;
            pars.epsilon = subtolerance; % Becca used 1e-5;
            if tau>0
                subsolution = denoise_bound(step,tau./alpha,-mu,Inf,pars);
            else
                subsolution = step.*(step>0);
            end
    end           
end

% =====================================
% = Termination Criteria Computation: =
% =====================================
function converged = checkconvergence(iter,miniter,stopcriterion,tolerance,...
                        dx, x, cputime, objective)
	converged = 0;
    if iter >= miniter %no need to check if miniter not yet exceeded

        switch stopcriterion
            case 1
                % Simply exhaust the maximum iteration budget
                converged = 0;
            case 2
                % Terminate after a specified CPU time (in seconds)
                converged = (cputime >= tolerance);
            case 3
                % Relative changes in iterate
                rr = ((sum(dx(:).^2)./sum(x(:).^2)))
                converged = ((sum(dx(:).^2)./sum(x(:).^2)) <= tolerance^2);
            case 4
                % relative changes in objective
                converged = (( abs( objective(iter+1) - objective(iter))...
                    ./abs(objective(iter)) ) <= tolerance);
            case 5
                % complementarity condition
                todo
            case 6
                % Norm of lagrangian gradient
                todo
        end
    end
end
 
function todo
    error('This function is not yet implemented, please be patient!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% End SPIRAL SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraints Section %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% those are used for spiral 

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

function [f_subs]= Projct_3D(f_vec)

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

function [ f_subs] = Alt2D_Proj(f_vec,maxit) 
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
              if   0 <= c_new && c_new <= 1 && 0 <= p_new && p_new <= 1 && 0<=gp_new && gp_new<=1 && c_new <= p_new && p_new<=gp_new  
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%% Plot functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ FPRv, Recall,precision] = NOVroc_gen_PPV( fhat_recon, f_true, thresh )

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
    Recall = TPRv; 
    
    %Positive Predictive Value: (TP)/(TP+FP)
    PPV(i) = T_p/(T_p+F_p);
    precision = PPV; 

end

end

function [Recall,Accur,precision] = metrics(f_true,fhat_recon,thresh)
f_thresh = fhat_recon >= thresh;
TP = sum((f_thresh == 1) & (f_true == 1)); 
FP = sum((f_thresh == 1) & (f_true == 0));
TN = sum((f_thresh == 0) & (f_true == 0));
FN = sum((f_thresh == 0) & (f_true == 1));
Accur = (TP+TN)./(TP+FP+TN+FN); 
Recall = TP/(TP + FN);
precision = TP/(TP + FP);
end



