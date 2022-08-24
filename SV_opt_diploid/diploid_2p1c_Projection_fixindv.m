function [ f_subs ] = diploid_2p1c_Projection_fixindv(f_vec,maxit)
% This is the function to alternate between the three 2D projections for the
% 2 parent 1 child diploid model (MeMeA 2020?)
% First we fix the parent's signals and let the child vary, then we look at
% two other 2D regions where we fix one of the parents and the child and
% let the remaining parent vary.


% determine length of vector n (number of observations)
N = length(f_vec);
n = N/6;

% unshuffle binary vectors z=(z_C; z_F; z_M), y=(y_C; y_F; y_M), where z_i
% are the homozygous signals and y_i are the heterozygous signals
zC_vec = f_vec(1:n);
yC_vec = f_vec((n+1):2*n);
zF_vec = f_vec((2*n+1):3*n);
yF_vec = f_vec((3*n+1):4*n);
zM_vec = f_vec((4*n+1):5*n);
yM_vec = f_vec((5*n+1):6*n);


%OLD
% z=f_vec(1:N/2);
% y=f_vec((N/2 + 1): N);
% 
% zC_vec = z(1:n);
% zF_vec = z((n+1):2*n);
% zM_vec = z((2*n+1):3*n);
% 
% yC_vec = y(1:n);
% yF_vec = y((n+1):2*n);
% yM_vec = y((2*n+1):3*n);
%

% Initialize projections
for i = 1:n
    
    zC=zC_vec(i);
    yC=yC_vec(i);
    zF=zF_vec(i);
    yF=yF_vec(i);
    zM=zM_vec(i);
    yM=yM_vec(i);
    
    % start counter j for subproblem iterations
    j = 1;
    %thresh = 0.01;
    
    while (j <= maxit)
        
%         if j==1
%             disp(i)
%         end
        
        %----------------------------------------------------------------------
        % Fix zF, yF and zM, yM for each location n.  Calculate zC and yC.
        %----------------------------------------------------------------------
        
        if (j == 1)
            
            % Check if parent-child observation already satisfies
            % constraints.
            
            if(0 <= zF && zF <= 1 && 0 <= zM && zM <= 1  && 0 <= zC && zC <= 1 && ...
                    0 <= yF && yF <= 1 && 0 <= yM && yM <= 1 && 0 <= yC && yC <= 1 && ...
                    0 <= (yF + zF) && (yF + zF) <= 1 && 0 <= (yM + zM) && (yM + zM) <= 1 && ...
                    0 <= (yC + zC) && (yC + zC) <= 1 && 0 <= zC  && zC <= (zF + yF) &&...
                    0 <= zC && zC <= (zM + yM) && 0 <= zF && zF <= (zC + yC) && ...
                    0 <= zM && zM <= (zC + yC) && max(zF + zM - 1, 0) <= zC &&...
                    yC <= min(zF + yF + zM + yM, 1))
                
                
                j = maxit;
            end
            
            %            % If constraints are not satisfied, then initialize projections
            
            %            zF=0.5;
            %            yF=0.5;
            %            zM=0.5;
            %            yM=0.5;
            
            % If constraints are not satisfied, then initialize projections
            % by thresholding,
            
            zFpos = sort([0,zF,1]);
            zF = zFpos(2);
            zMpos = sort([0,zM,1]);
            zM = zMpos(2);
            
            yFpos = sort([0,yF,1]);
            yF = yFpos(2);
            yMpos = sort([0,yM,1]);
            yM = yMpos(2);
            
            % Enforce constraint that individual het. + hom. SV <= 1.
            if ( (zF+yF) > 1 )
                zF = 0.5;
                yF = 0.5;
            end
            if ( (zM+yM) > 1 )
                zM = 0.5;
                yM = 0.5;
            end
            
            
        else
            % Continue projections w/ previous iterate.
            zF = zF_new;
            yF = yF_new;
            zM = zM_new;
            yM = yM_new;
            
        end
        
        % Begin projections to calculate zC and yC
        %disp(j)
        [zC_new, yC_new] = Diploid_2p1c_FS1_REVISED(zC, yC, zF,yF,zM,yM);
        
        %          disp(zC_new)
        %          disp(yC_new)
        
        %----------------------------------------------------------------------
        % Given updated zC and yC, fix zF, yF, zC, yC and calculate new zM and
        % yM for same location.
        %----------------------------------------------------------------------
        
        zC = zC_new;
        yC = yC_new;
        
        % Enter the child signals, followed by the fixed parents signals,
        % followed by the parent's signals that need to be projected
        [zM_new, yM_new] = Diploid_2p1c_FS2_clean(zC,yC,zF,yF,zM,yM);
        %          disp("M")
        
        %----------------------------------------------------------------------
        % Given updated zM and yM, fix zM, yM, zC, yC and calculate new zF and
        % yF for same location.
        %----------------------------------------------------------------------
        
        zM = zM_new;
        yM = yM_new;
        
        % Enter the child signals, followed by the fixed parents signals,
        % followed by the parent's signals that need to be projected
        [zF_new, yF_new] = Diploid_2p1c_FS2_clean(zC,yC,zM,yM,zF,yF);
        %          disp("F")
        
        % Increase counter
        j = j+1;
        
    end
    
    
    zC_vec(i)=zC_new;
    yC_vec(i)=yC_new;
    zF_vec(i)=zF_new;
    yF_vec(i)=yF_new;
    zM_vec(i)=zM_new;
    yM_vec(i)=yM_new;
    
    
end

f_subs = [zC_vec;yC_vec;zF_vec;yF_vec;zM_vec;yM_vec];
%f_subs = [zC_vec;zF_vec;zM_vec;yC_vec;yF_vec;yM_vec];
