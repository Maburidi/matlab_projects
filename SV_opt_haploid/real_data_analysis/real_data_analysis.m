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

f_c = True_sig{:,11};
f_c(f_c>0)=1;
f_p = True_sig{:,12};
f_p(f_p>0)=1;
f_gp = True_sig{:,8};
f_gp(f_gp>0)=1;

f_p2 = True_sig{:,14};
f_p2(f_p2>0)=1;
f_gp2 = True_sig{:,1};
f_gp2(f_gp2>0)=1;



y_c = Observs{:,11};
y_p = Observs{:,12};
y_gp = Observs{:,8};
%%
figure(1)

subplot(1,3,1)
hist(y_gp(find(f_gp ==1)),0:100)
xlim([0 15])
ylim([0 70])
xticks(0:15)

title("Histogram of the Corresponding $y_{gp}$ to SVs in $f_{gp}$",Interpreter="latex",FontSize=15)
xlabel("Number of Fregments")
ylabel("Counts")

subplot(1,3,2)
hist(y_p(find(f_p ==1)),0:100)
xlim([0 15])
ylim([0 70])
xticks(0:15)
title("Histogram of the Corresponding $y_{p}$ to SVs in $f_p$",Interpreter="latex",FontSize=15)
xlabel("Number of Fregments")
ylabel("Counts")

subplot(1,3,3)
hist(y_c(find(f_c ==1)),0:100)
xlim([0 15])
ylim([0 70])
xticks(0:15)

title("Histogram of the Corresponding $y_{c}$ to SVs in $f_c$",Interpreter="latex",FontSize=15)
xlabel("Number of Fregments")
ylabel("Counts")

%%

ind_c= find(f_c==1);
ind_p= find(f_p==1);
ind_gp= find(f_gp==1);
ind_ = [ind_c ;ind_p; ind_gp]; 
ind = unique(ind_); 
m=0;
n=0;
ind_vio=[]; 
for i=1:length(ind)
 if f_c(ind(i)) <= f_p(ind(i)) && f_p(ind(i)) <= f_gp(ind(i))  
     m=m+1;  
 else 
     n=n+1; 
     ind_vio = [ind_vio ind(i)];
 end 
end 

figure(2)
hist([ones(m,1); 2*ones(n,1)],1:3)
xticklabels(["Satisfy" , "Not Satisfy"])
title("number of SVs satisfies the constrains",Interpreter="latex",FontSize=10)
%xlabel("Number of SVs")
ylabel("Count")

%%

mm=0;
nn=0;
zz=0;
z=0; 
m=0;
n=0; 
b=0; 
ind_vio_c1 = [];
ind_vio_c2 = [];
ind_vio_c3 = [];

for i=1:length(ind)
 if f_c(ind(i)) <= f_p(ind(i)) 
     mm=mm+1;  
 end

 if f_p(ind(i)) <= f_gp(ind(i))
     nn=nn+1; 
 end 

 if  f_c(ind(i)) <= f_gp(ind(i)) 
    zz=zz+1;
 end 
 
 if f_c(ind(i)) > f_p(ind(i)) 
     m=m+1;
     ind_vio_c1 = [ind_vio_c1 ind(i)];

 end 

 if f_p(ind(i)) > f_gp(ind(i)) 
     n=n+1;
     ind_vio_c2 = [ind_vio_c2 ind(i)];

 end 

 if f_c(ind(i)) > f_gp(ind(i)) 
     z=z+1;
     ind_vio_c3 = [ind_vio_c3 ind(i)];
 end 

end 

figure(3)
%hist([ones(m,1); 2*ones(n,1)],1:3)
y = [mm m;nn n;zz z]
bar(y)
xticklabels(["f_c <= f_p" , "f_p <= f_{gp}", "f_c <= f_{gp}"])
title("number of SVs satisfies the constrains",Interpreter="latex",FontSize=15)
%xlabel("Number of SVs")
ylabel("Count")
legend("Satisfy" , "Not Satisfy")

%%
ind_p2 = f_p2(ind_vio);
ind_gp2 = f_gp2(ind_vio);
figure(4)
yy= [length(find(ind_p2==1))  length(find(ind_p2==0));length(find(ind_gp2==1)) length(find(ind_gp2==0)) ];
bar(yy)
xticklabels(["Parent" , "Grandparent"])
%title("number of SVs satisfies the constrains",Interpreter="latex",FontSize=15)
%xlabel("Number of SVs")
ylabel("Count")
legend("Ones" , "Zeros")


%---------------------

ind_p2_1 = f_p2(ind_vio_c1);
ind_p2_2 = f_p2(ind_vio_c2);

figure(5)
yy= [length(find(ind_p2_1==1))  length(find(ind_p2_1==0));length(find(ind_p2_2==1)) length(find(ind_p2_2==0)) ];
bar(yy)
xticklabels(["f_c <= f_p" , "f_p <= f_{gp}"])
%title("number of SVs satisfies the constrains",Interpreter="latex",FontSize=15)
%xlabel("Number of SVs")
ylabel("Count")
legend("Ones" , "Zeros")

%---------------------

ind_gp2_1 = f_gp2(ind_vio_c2);
ind_gp2_2 = f_gp2(ind_vio_c3);
figure(6)
yy= [length(find(ind_gp2_1==1))  length(find(ind_gp2_1==0));length(find(ind_gp2_2==1)) length(find(ind_gp2_2==0)) ];
bar(yy)
xticklabels(["f_p <= f_{gp}" , "f_c <= f_{gp}"])
%title("number of SVs satisfies the constrains",Interpreter="latex",FontSize=15)
%xlabel("Number of SVs")
ylabel("Count")
legend("Ones" , "Zeros")














