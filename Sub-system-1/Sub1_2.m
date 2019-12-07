clc;clear;
close all;
tic
%%
Material_Liner = readmatrix('material_liner.csv','Range','B2:G38');
Material_Shell = readmatrix('material_shell.csv','Range','B2:G17');
%%
[~,idx] = sort(Material_Liner(:,2)); 
Material_Liner = Material_Liner(idx,:);
[~,idx] = sort(Material_Shell(:,2)); 
Material_Shell = Material_Shell(idx,:);

% rho_s = Material_Shell(:,2);
% ESE_s= Material_Shell(:,3);
% EE = Material_Liner(:,5);
% CO2e = Material_Liner(:,6);
% EE_s= Material_Shell(:,5);
% CO2e_s = Material_Shell(:,6);
% figure(1)
% plot(rho_s,EE_s,'o')
% figure(2)
% plot(rho_s,CO2e_s,'o')
% figure(3)
% plot(rho_s,ESE_s,'o')

TF1 = Material_Liner(:,5)>1.5e8;
TF2 = Material_Liner(:,6)>10;
TFall1 = TF1 & TF2;
Material_Liner(TFall1,:)=[];

TF3 = Material_Shell(:,5)>1.2e8;
TF4 = Material_Shell(:,6)>5;
TFall2 = TF3 & TF4;
Material_Shell(TFall2,:)=[];

rho_l = Material_Liner(:,2);
sigma_l = Material_Liner(:,3);
EE_l = Material_Liner(:,5);
CO2e_l = Material_Liner(:,6);

global Beta1 Beta2 Beta3
Beta1 = polyfit(rho_l,sigma_l,2);
Beta2 = polyfit(rho_l,EE_l,1);
Beta3 = polyfit(rho_l,CO2e_l,1);

rho_s = Material_Shell(:,2);
ESE_s = Material_Shell(:,3);
EE_s = Material_Shell(:,5);
CO2e_s = Material_Shell(:,6);

global Beta4 Beta5 Beta6
Beta4 = polyfit(rho_s,ESE_s,2);
Beta5 = polyfit(rho_s,EE_s,1);
Beta6 = polyfit(rho_s,CO2e_s,1);
%%
global v0 mh r
v0 = 5;
mh = 5;
r = 0.08;
%%
rng default
% PLA = (v0./M.^(0.5)).*(2.*pi.*R.*sigma).^0.5;
% dhmax = v0.*(M./(2.*pi.*R.*sigma)).^0.5;
x0 = [0.02 22 0.005 1050];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.02 22 0.005 1050];
ub = [0.02 860 0.01 1580];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','Plotfcn','optimplotfval');
[xmin,fval] = fmincon(@PLA,x0,A,b,Aeq,beq,lb,ub,@dhmax,options);
toc
%%
function pla = PLA(x)
global Beta1 Beta4 v0 mh r
R1 = r+x(1);
R2 = r+x(1)+x(3);
sigma = Beta1(1).*x(2).^2+Beta1(2).*x(2)+Beta1(3);
ESE = Beta4(1).*x(4).^2+Beta4(2)+Beta4(3);
M = mh+2/3*pi.*(x(2).*(R1.^3-r^3)+x(4).*(R2.^3-R1.^3));
v1 = (v0^2-(4.*ESE.*pi^3.*R2.^2.*x(3).^2.*x(3))./M).^0.5;
pla = v1./M.^0.5.*(2.*pi.*(R1).*(sigma)).^0.5/9.81;
end

function [c,ceq] = dhmax(x)
global Beta1 Beta2 Beta3 Beta4 Beta5 Beta6 v0 mh r
R1 = r+x(1);
R2 = r+x(1)+x(3);
sigma = Beta1(1).*x(2).^2+Beta1(2).*x(2)+Beta1(3);
ESE = Beta4(1).*x(4).^2+Beta4(2)+Beta4(3);
EE_l = Beta2(1).*x(2)+Beta2(2);
CO2e_l = Beta3(1).*x(2)+Beta3(2);
EE_s = Beta5(1).*x(4)+Beta5(2);
CO2e_s = Beta6(1).*x(4)+Beta6(2);
M_l = 2/3*pi.*x(2).*(R1.^3-r^3);
M_s =  2/3*pi.*x(4).*(R2.^3-R1.^3);
M = mh+M_l+M_s;
v1 = (v0^2-(4.*ESE.*pi^3.*R2.^2.*x(3).^2.*x(3))./M).^0.5;
c(1) = v1.*(M./(2.*pi.*R1.*sigma)).^0.5 - 0.05;
c(2) = M_l+M_s-0.5;
c(3) = M_l.*EE_l+M_s.*EE_s-15e7; 
c(4) = M_l.*CO2e_l+M_s.*CO2e_s-10;
ceq = [];
end