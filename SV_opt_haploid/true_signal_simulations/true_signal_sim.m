% This sript detect the SV on GP-> P-> C  Diploid case
% =====================================================================    
% =    True SV vectors and obsevration vectors Simulations    =
% =====================================================================
% this piece of code is to generate the vector f for GP, f for P and f for C
clc
clear
close all 

n = 10^5;        % length of the each signal (f vector)

%---- k_gp1 >= k_gp2
k_gp1 = round(0.05*n);          % number of SV in GP1, 5% of the signal
k_gp2 = round(0.05*n);            % number of SV in GP2, 5% of the signal

%k_p = round(0.5*k_gp);    % number of SV in P
%k_c = round(0.25*k_gp);   % number of SV in C

p = 0.0;              % percent of II copy SVs in gp1 and gp2, if p1 = 1, then the parent is completely homozygous, two copies 
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


