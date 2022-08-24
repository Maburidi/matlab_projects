function [f_subs] = consts_diploid_1gp1p1c(f_vec , maxit)

% determine length of vector n (number of observations)
N = length(f_vec);
n = N/6;

% unshuffle binary vectors z=(z_C; z_F; z_M), y=(y_C; y_F; y_M), where z_i
% are the homozygous signals and y_i are the heterozygous signals
zgp_vec = f_vec(1:n);
ygp_vec = f_vec((n+1):2*n);
zp_vec = f_vec((2*n+1):3*n);
yp_vec = f_vec((3*n+1):4*n);
zc_vec = f_vec((4*n+1):5*n);
yc_vec = f_vec((5*n+1):6*n);


for i = 1:n
    
    zC=zc_vec(i);
    yC=yc_vec(i);
    zP=zp_vec(i);
    yP=yp_vec(i);
    zGP=zgp_vec(i);
    yGP=ygp_vec(i);
    
    % start counter j for subproblem iterations
    j = 1;
    
    while (j < maxit)
         % Check if parent-child observation already satisfies constraints.

         if(0 <= zC && zC <= 1 && 0 <= zP && zP <= 1  && 0 <= zGP && zGP <= 1 && ...
                    0 <= yP && yP <= 1 && 0 <= yGP && yGP <= 1 && 0 <= yC && yC <= 1 && ...
                    0 <= (yP + zP) && (yP + zP) <= 1 && ...
                    0 <= (yGP + zGP) && (yGP + zGP) <= 1 && ...
                    0 <= (yC + zC) && (yC + zC) <= 1 && ...
                    zC <= (zP + yP) &&...
                    zP <= (zGP + yGP) && ...
                    zP <= (zC + yC) &&...
                    zGP <= (zP + yP))
                 
                    zC_new = zC;
                    yC_new = yC;
                    zP_new = zP;
                    yP_new = yP;
                    zGP_new = zGP;
                    yGP_new = yGP; 

                j = maxit;
                continue 
        else 

        if j==1            %%%%%  initialization 
            zPpos = sort([0,zP,1]);
            zP = zPpos(2);
            zGPpos = sort([0,zGP,1]);
            zGP = zGPpos(2);
            
            yPpos = sort([0,yP,1]);
            yP = yPpos(2);
            yGPpos = sort([0,yGP,1]);
            yGP = yGPpos(2);
            % Enforce constraint that individual het. + hom. SV <= 1.
            if ((zP+yP)>1)
                zP = 0.5;
                yP = 0.5;
            end
            if ((zGP+yGP)>1)
                zGP = 0.5;
                yGP = 0.5;
            end
        end  
        end 

        % Begin projections 
        % Project zC and yC 
        [zC_new, yC_new] = Diploid_1gp1p1c_CASE1(zC, yC, zP,yP,zGP,yGP);
        zC = zC_new;
        yC = yC_new;
        
        % Project zP and yP 
        [zP_new, yP_new] = Diploid_1gp1p1c_ParentCase(zP,yP,zGP,yGP,zC,yC); 
             zP = zP_new;
             yP = yP_new;
%         if zGP >= zC
%              [zP_new, yP_new] = Diploid_1gp1p1c_CASE1(zP,yP,zGP,yGP,zC,yC);
%              zP = zP_new;
%              yP = yP_new;
%         else 
%              [zP_new, yP_new] = Diploid_1gp1p1c_CASE2(zP,yP,zGP,yGP,zC,yC);         
%              zP = zP_new;
%              yP = yP_new;
%         end  
       
        % Project zGP and yGP 
        [zGP_new, yGP_new] = Diploid_1gp1p1c_CASE1(zGP,yGP,zP,yP,zC,yC);
        zGP = zGP_new; 
        yGP = yGP_new; 
         

        j = j+1;
    end 
  
    zc_vec(i)=zC_new;
    yc_vec(i)=yC_new;
    zp_vec(i)=zP_new;
    yp_vec(i)=yP_new;
    zgp_vec(i)=zGP_new;
    ygp_vec(i)=yGP_new;

    
end 
f_subs = [zgp_vec;ygp_vec;zp_vec;yp_vec;zc_vec;yc_vec];

end 

function [zC_new, yC_new]= Diploid_1gp1p1c_CASE1(zC,yC,zP,yP,zGP,yGP)



% This is the projections for the child's signals (zC,yC) if we fix both
% parents signals (zF,zM,yF,yM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0 <= zP && zP <= 1 && 0 <= zGP && zGP <= 1 && ...
            0 <= yP && yP <= 1 && 0 <= yGP && yGP <= 1 && ...
            0 <= (yP + zP) && (yP + zP) <= 1 &&...
            0 <= (yGP + zGP) && (yGP + zGP) <= 1)
        %fits constraits so we are ok
else
    disp("ERROR in FS1 Projections -- Not feasible")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if yP ==1 && zP ==0    %% 1- project to triangle
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif yC >= 1-zC && yC >= zC-1 && yC <= zC+1
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 3")
                    zC_new=zC;
                    yC_new=yC;
                end

elseif yP ==0 && zP ==0   %% 2- project to a line 
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= yC && yC <= 1 && zC >= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC >= 1 
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 3")
                    zC_new=zC;
                    yC_new=yC;
                end

elseif 0 < yP && yP < 1 && zP ==0    %% 3- 
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= yP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= yP
                    zC_new=yP;
                    yC_new=0;
                elseif zC >= yP && 0 <= yC && yC <= 1-yP
                    zC_new=yP;
                    yC_new=yC;
                elseif 1-yP <= yC && yC <= zC+1-2*yP
                    zC_new=yP;
                    yC_new=1-yP;
                elseif zC+1-2*yP <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
%                     %disp("Interior 2")
                    zC_new=zC;
                    yC_new=yC;
                end

elseif  yP == 0 && zP ==1       %% 4- 

               if yC >= zC+1
                      zC_new=0;
                      yC_new=1;
                elseif zC-1 <= yC && yC <= zC+1 && yC <= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                elseif yC <= zC-1
                      zC_new=0;
                      yC_new=1;
                elseif zC-1 <= yC && yC <= zC+1 && yC >= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                else
                      %disp("Interior 24")
                      zC_new=zC;
                      yC_new=yC;
                end

%[5]
elseif  0 < yP && yP < 1 && 0 < zP && zP < 1 && zP+yP ==1     %% 5- 

                if zP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zP <= yC && yC <= zP
                    zC_new=0;
                    yC_new=zP;
                elseif zC - zP <= yC && yC <= zC + zP && yC <= zP - zC
                    zC_new=(1/2)*(zP-yC+zC);
                    yC_new=(1/2)*(zP+yC-zC);
                elseif zC <= zP && yC <= zC - zP
                    zC_new=zP;
                    yC_new=0;
                elseif zP <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif zC-1 <= yC && yC <= zC+1 && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end

%[6]
elseif  0 < yP && yP < 1 && 0 < zP && zP < 1 && zP+yP <1     
                if zP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zP <= yC && yC <= zP
                    zC_new=0;
                    yC_new=zP;
                elseif zC - zP <= yC && yC <= zC + zP && yC <= zP - zC
                    zC_new=(1/2)*(zP-yC+zC);
                    yC_new=(1/2)*(zP+yC-zC);
                elseif zC <= zP && yC <= zC - zP
                    zC_new=zP;
                    yC_new=0;
                elseif zP <= zC && zC <= zP+yP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= zP+yP
                    zC_new=zP+yP;
                    yC_new=0;
                elseif zC >= yP+zP && 0 <= yC && yC <= 1-(yP+zP)
                    zC_new=yP+zP;
                    yC_new=yC; 
                elseif 1-(yP+zP) <= yC && yC <= zC+1-2*(yP+zP)
                    zC_new=yP+zP;          
                    yC_new=1-(yP+zP);
                elseif zC+1-2*(yP+zP) <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);         
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end
 
%[7]
elseif  yP == 0 && 0 < zP && zP < 1        %% 7-   
                if zP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zP <= yC && yC <= zP
                    zC_new=0;
                    yC_new=zP;
                elseif zC - zP <= yC && yC <= zC + zP && yC <= zP - zC
                    zC_new=(1/2)*(zP-yC+zC);
                    yC_new=(1/2)*(zP+yC-zC);
                elseif zC <= zP && yC <= zC - zP
                    zC_new=zP;
                    yC_new=0;
                elseif yC <= 0 && zC >= zP
                    zC_new=zP;
                    yC_new=0;
                elseif zC >= zP && 0 <= yC && yC <= 1-zP
                    zC_new=zP;
                    yC_new=yC; 
                elseif 1-zP <= yC && yC <= zC+1-2*zP
                    zC_new= zP;          
                    yC_new= 1-zP;
                elseif zC+1-2*zP <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);         
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end
    
end 

end 


function [zC_new, yC_new]= Diploid_1gp1p1c_CASE2(zC,yC,zP,yP,zGP,yGP)

% This is the projections for the child's signals (zC,yC) if we fix both
% parents signals (zF,zM,yF,yM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if (0 <= zP && zP <= 1 && 0 <= zGP && zGP <= 1 && ...
             0 <= yP && yP <= 1 && 0 <= yGP && yGP <= 1 && ...
             0 <= (yP + zP) && (yP + zP) <= 1 &&...
             0 <= (yGP + zGP) && (yGP + zGP) <= 1)
         %fits constraits so we are ok
 else
     disp("ERROR in FS1 Projections -- Not feasible")
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [1] zGP = 0 , zP =0, yP =0 
if zGP == 0 && zP == 0 && yP == 0          %% 2- project to a line 
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= yC && yC <= 1 && zC >= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC >= 1 
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 3")
                    zC_new=zC;
                    yC_new=yC;
                end
%[2]
elseif zGP == 0 && zP == 0 && yP == 1    %% project to triangle
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif yC >= 1-zC && yC >= zC-1 && yC <= zC+1
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 3")
                    zC_new=zC;
                    yC_new=yC;
                end

 %[3]               
elseif zGP == 0 && zP == 1 && yP == 0    %% project to triangle
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif yC >= 1-zC && yC >= zC-1 && yC <= zC+1
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 3")
                    zC_new=zC;
                    yC_new=yC;
                end

 %[4]                
elseif zGP == 1 && zP == 0 && yP == 0    %% project to dot 
                    zC_new=0;
                    yC_new=1;

%[5]
elseif  zGP == 1 && zP == 0 && yP == 1    % line 
               if yC >= zC+1
                      zC_new=0;
                      yC_new=1;
                elseif zC-1 <= yC && yC <= zC+1 && yC <= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                elseif yC <= zC-1
                      zC_new=0;
                      yC_new=1;
                elseif zC-1 <= yC && yC <= zC+1 && yC >= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                else
                      %disp("Interior 24")
                      zC_new=zC;
                      yC_new=yC;
                end

%[6]
elseif  zGP == 1 && zP == 1 && yP == 0      % line 
               if yC >= zC+1
                      zC_new=0;
                      yC_new=1;
                elseif zC-1 <= yC && yC <= zC+1 && yC <= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                elseif yC <= zC-1
                      zC_new=0;
                      yC_new=1;
                elseif zC-1 <= yC && yC <= zC+1 && yC >= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                else
                      %disp("Interior 24")
                      zC_new=zC;
                      yC_new=yC;
               end
%----------------------------------
%[7]
elseif zGP == 0 && zP ==0  && 0 < yP && yP < 1   
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= yP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= yP
                    zC_new=yP;
                    yC_new=0;
                elseif zC >= yP && 0 <= yC && yC <= 1-yP
                    zC_new=yP;
                    yC_new=yC;
                elseif 1-yP <= yC && yC <= zC+1-2*yP
                    zC_new=yP;
                    yC_new=1-yP;
                elseif zC+1-2*yP <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
%                     %disp("Interior 2")
                    zC_new=zC;
                    yC_new=yC;
                end
               
%[8]
elseif zGP == 0 && 0 < zP && zP < 1  && yP ==0    
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= zP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= zP
                    zC_new=zP;
                    yC_new=0;
                elseif zC >= zP && 0 <= yC && yC <= 1-zP
                    zC_new=zP;
                    yC_new=yC;
                elseif 1-zP <= yC && yC <= zC+1-2*zP
                    zC_new=zP;
                    yC_new=1-zP;
                elseif zC+1-2*zP <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
%                     %disp("Interior 2")
                    zC_new=zC;
                    yC_new=yC;
                end

%[9]
elseif zGP == 0 && 0 < zP && zP < 1  && 0 < yP && yP < 1 && yP + zP <1     
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= zP+yP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= zP+yP
                    zC_new=zP+yP;
                    yC_new=0;
                elseif zC >= zP+yP && 0 <= yC && yC <= 1-(zP+yP)
                    zC_new=zP+yP;
                    yC_new=yC;
                elseif 1-(zP+yP) <= yC && yC <= zC+1-2*(zP+yP)
                    zC_new=zP+yP;
                    yC_new=1-(zP+yP);
                elseif zC+1-2*(zP+yP) <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
%                     %disp("Interior 2")
                    zC_new=zC;
                    yC_new=yC;
                end

%[10]
elseif zGP == 0 && 0 < zP && zP < 1  && 0 < yP && yP < 1 && yP + zP ==1     
                if 0 <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif yC <= 0 && zC <= 0
                    zC_new=0;
                    yC_new=0;
                elseif 0 <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif yC >= 1-zC && yC >= zC-1 && yC <= zC+1
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >=zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 3")
                    zC_new=zC;
                    yC_new=yC;
                end
%[11]
elseif  zGP == 1 && zP == 0  && 0 < yP && yP < 1 

               if yC >= zC+1
                      zC_new=0;
                      yC_new=1;
                elseif zC+1 - 2*yP <= yC && yC <= zC+1 && yC <= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                elseif yC <= zC+1 - 2*yP
                      zC_new=yP;
                      yC_new=1-yP;
                elseif (zC+1 - 2*yP) <= yC && yC <= zC+1 && yC >= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                else
                      %disp("Interior 24")
                      zC_new=zC;
                      yC_new=yC;
                end


%[12]
elseif  zGP == 1 && 0 < zP && zP < 1 && yP == 0
               if yC >= zC+1
                      zC_new=0;
                      yC_new=1;
                elseif zC+1 - 2*zP <= yC && yC <= zC+1 && yC <= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                elseif yC <= zC+1 - 2*zP
                      zC_new=zP;
                      yC_new=1-zP;
                elseif (zC+1 - 2*zP) <= yC && yC <= zC+1 && yC >= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                else
                      %disp("Interior 24")
                      zC_new=zC;
                      yC_new=yC;
                end

%[13]
elseif  zGP == 1 && 0 < zP && zP < 1 && 0 < yP && yP < 1
               if yC >= zC+1
                      zC_new=0;
                      yC_new=1;
                elseif zC+1 - 2*(zP+yP) <= yC && yC <= zC+1 && yC <= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                elseif yC <= zC+1 - 2*(zP+yP)
                      zC_new=zP+yP;
                      yC_new=1-(zP+yP);
                elseif (zC+1 - 2*(zP+yP)) <= yC && yC <= zC+1 && yC >= 1-zC
                      zC_new=(1/2)*(1-yC+zC);
                      yC_new=(1/2)*(1+yC-zC);
                else
                      %disp("Interior 24")
                      zC_new=zC;
                      yC_new=yC;
                end


%[14]
elseif  0 < zGP && zGP < 1  &&  zP==0 && yP == 0
    if yC >=1 
         zC_new=0;
         yC_new=1;
    elseif yC <= zGP
         zC_new=0;
         yC_new=zGP;
    else  
         zC_new=0;
         yC_new=yC;
    end 

%[15]
elseif  0 < zGP && zGP < 1  &&  zP==0 && yP == 1
                if zGP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zGP <= yC && yC <= zGP
                    zC_new=0;
                    yC_new=zGP;
                elseif zC - zGP <= yC && yC <= zC + zGP && yC <= zGP - zC
                    zC_new=(1/2)*(zGP-yC+zC);
                    yC_new=(1/2)*(zGP+yC-zC);
                elseif zC <= zGP && yC <= zC - zGP
                    zC_new=zGP;
                    yC_new=0;
                elseif zGP <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif zC-1 <= yC && yC <= zC+1 && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end

%[16]
elseif  0 < zGP && zGP < 1  &&  zP==1 && yP == 0
                if zGP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zGP <= yC && yC <= zGP
                    zC_new=0;
                    yC_new=zGP;
                elseif zC - zGP <= yC && yC <= zC + zGP && yC <= zGP - zC
                    zC_new=(1/2)*(zGP-yC+zC);
                    yC_new=(1/2)*(zGP+yC-zC);
                elseif zC <= zGP && yC <= zC - zGP
                    zC_new=zGP;
                    yC_new=0;
                elseif zGP <= zC && zC <= 1 && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= zC-1 && zC >= 1
                    zC_new=1;
                    yC_new=0;
                elseif zC-1 <= yC && yC <= zC+1 && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end

% [17]             
elseif  0 < zGP && zGP < 1  &&  zP==0 && 0 < yP && yP < 1
                if zGP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zGP <= yC && yC <= zGP
                    zC_new=0;
                    yC_new=zGP;
                elseif zC - zGP <= yC && yC <= zC + zGP && yC <= zGP - zC
                    zC_new=(1/2)*(zGP-yC+zC);
                    yC_new=(1/2)*(zGP+yC-zC);
                elseif zC <= zGP && yC <= zC - zGP
                    zC_new=zGP;
                    yC_new=0;
                elseif zGP <= zC && zC <= yP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= yP
                    zC_new=yP;
                    yC_new=0;
                elseif zC >= yP && 0 <= yC && yC <= 1-yP
                    zC_new=yP;
                    yC_new=yC;
                elseif 1-yP <= yC && yC <= zC+1-2*yP
                    zC_new=yP;          
                    yC_new=1-yP;
                elseif zC+1-2*yP <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);  
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end

%[18]
elseif  0 < zGP && zGP < 1  &&  0 < zP && zP < 1 &&  0 < yP && yP < 1
                if zGP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zGP <= yC && yC <= zGP
                    zC_new=0;
                    yC_new=zGP;
                elseif zC - zGP <= yC && yC <= zC + zGP && yC <= zGP - zC
                    zC_new=(1/2)*(zGP-yC+zC);
                    yC_new=(1/2)*(zGP+yC-zC);
                elseif zC <= zGP && yC <= zC - zGP
                    zC_new=zGP;
                    yC_new=0;
                elseif zGP <= zC && zC <= zP && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= zP
                    zC_new=zP;
                    yC_new=0;
                elseif zC >= zP && 0 <= yC && yC <= 1-zP
                    zC_new=zP;
                    yC_new=yC;
                elseif 1-zP <= yC && yC <= zC+1-2*zP
                    zC_new=zP;          
                    yC_new=1-zP;
                elseif zC+1-2*zP <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);  
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end
                     
%[19]
elseif  0 < zGP && zGP < 1  &&  0 < zP && zP < 1 &&  yP == 0
                if zGP <= yC && yC <= 1 && zC <= 0
                    zC_new=0;
                    yC_new=yC;
                elseif zC + zGP <= yC && yC <= zGP
                    zC_new=0;
                    yC_new=zGP;
                elseif zC - zGP <= yC && yC <= zC + zGP && yC <= zGP - zC
                    zC_new=(1/2)*(zGP-yC+zC);
                    yC_new=(1/2)*(zGP+yC-zC);
                elseif zC <= zGP && yC <= zC - zGP
                    zC_new=zGP;
                    yC_new=0;
                elseif zGP <= zC && zC <= (zP+yP) && yC <= 0
                    zC_new=zC;
                    yC_new=0;
                elseif yC <= 0 && zC >= zP+yP
                    zC_new=zP+yP;
                    yC_new=0;
                elseif zC >= zP+yP && 0 <= yC && yC <= 1-zP-yP
                    zC_new=zP+yP;
                    yC_new=yC;
                elseif 1-zP-yP <= yC && yC <= zC+1-2*(zP+yP)
                    zC_new=zP+yP;          
                    yC_new=1-zP-yP;
                elseif zC+1-2*(zP+yP) <= yC && yC <= 1+zC && yC >= 1-zC
                    zC_new=(1/2)*(1-yC+zC);
                    yC_new=(1/2)*(1+yC-zC);  
                elseif yC >= zC+1 && yC >=1
                    zC_new=0;
                    yC_new=1;
                else
                    %disp("Interior 11")
                    zC_new=zC;
                    yC_new=yC;
                end

else 
    disp("Projections are Not Complete")


               
end 

end 

 
function [zP_new, yP_new]= Diploid_1gp1p1c_ParentCase(zP,yP,zGP,yGP,zC,yC)

if zC + yC <= zGP + yGP && zGP <= zC
  [zP_new, yP_new] =  Diploid_1gp1p1c_CASE1(zP,yP,zC,yC,zGP,yGP); 

elseif zC + yC >= zGP + yGP && zGP >= zC
  [zP_new, yP_new] =  Diploid_1gp1p1c_CASE1(zP,yP,zGP,yGP,zC,yC); 

elseif zC + yC > zGP + yGP && zGP < zC
  [zP_new, yP_new] =  Diploid_1gp1p1c_CASE2(zP,yP,zC,yC,zGP,yGP); 

elseif zC + yC < zGP + yGP && zGP > zC
  [zP_new, yP_new] =  Diploid_1gp1p1c_CASE2(zP,yP,zGP,yGP,zC,yC); 

else 
    disp("Coding Error")

end 
end 


