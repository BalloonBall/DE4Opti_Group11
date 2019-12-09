clc;clear;
close all;
%% Data import
% import dataset from csv
Material_Liner = readmatrix('material_liner.csv','Range','B2:G38');
Material_Shell = readmatrix('material_shell.csv','Range','B2:G17');
%% Data processing
% sort dataset by density in asending order
[~,idx] = sort(Material_Liner(:,2)); 
Material_Liner = Material_Liner(idx,:);
[~,idx] = sort(Material_Shell(:,2)); 
Material_Shell = Material_Shell(idx,:);

% exclude materials with significant outlier properties(embodied energy &
% CO2 footprint)
TF1 = Material_Liner(:,5)>1.5e8;
TF2 = Material_Liner(:,6)>10;
TFall1 = TF1 & TF2;
Material_Liner(TFall1,:)=[];

% import material properties of liner and shell into arrays
% Liner
rho_l = Material_Liner(:,2);
sigma_l = Material_Liner(:,3);
EE_l = Material_Liner(:,5);
CO2e_l = Material_Liner(:,6);
% Shell
rho_s = Material_Shell(:,2);
ESE_s = Material_Shell(:,3);
EE_s = Material_Shell(:,5);
CO2e_s = Material_Shell(:,6);

% Variable linking:
% link yeild stress, embodied energy and CO2 footprint with density
global Beta1 Beta2 Beta3 Beta4 Beta5 Beta6
% Liner
Beta1 = polyfit(rho_l,sigma_l,2);
Beta2 = polyfit(rho_l,EE_l,1);
Beta3 = polyfit(rho_l,CO2e_l,1);
% Shell
Beta4 = polyfit(rho_s,ESE_s,2);
Beta5 = polyfit(rho_s,EE_s,2);
Beta6 = polyfit(rho_s,CO2e_s,2);
%% Define parameters
global v0 mh r
v0 = 5;
mh = 5;
r = 0.08;
%% Optimisation - Interior-point
tic
x0 = [0.02 22 0.005 1050];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.02 22 0.005 1050];
ub = [0.05 860 0.01 1580];
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','Plotfcn','optimplotfval');
[xmin_ip,fval_ip] = fmincon(@PLA,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
toc
%% Optimisation - SQP
tic
x0 = [0.02 22 0.005 1050];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.02 22 0.005 1050];
ub = [0.05 860 0.01 1580];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','Plotfcn','optimplotfval');
[xmin_sqp,fval_sqp] = fmincon(@PLA,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
toc
%% Optimisation - Active-set
tic
x0 = [0.02 22 0.005 1050];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.02 22 0.005 1050];
ub = [0.05 860 0.01 1580];
options = optimoptions('fmincon','Display','iter','Algorithm','active-set','Plotfcn','optimplotfval');
[xmin_as,fval_as] = fmincon(@PLA,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
toc