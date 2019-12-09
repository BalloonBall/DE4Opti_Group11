%% Nonlinear constraints
% c(1)maximum liner deflection
% c(2)total mass
% c(3)total embodied energy
% c(4)total cost of material

function [c,ceq] = nonlcon(x)
global Beta1 Beta2 Beta3 Beta4 Beta5 Beta6 v0 mh r
% Variable linking
R1 = r+x(1);
R2 = r+x(1)+x(3);
sigma = Beta1(1).*x(2).^2+Beta1(2).*x(2)+Beta1(3);
ESE = Beta4(1).*x(4).^2+Beta4(2).*x(4)+Beta4(3);
EE_l = Beta2(1).*x(2)+Beta2(2);
C_l = Beta3(1).*EE_l+Beta3(2);
EE_s = Beta5(1).*x(4)^2+Beta5(2).*x(4)+Beta5(3);
C_s = Beta6(1).*EE_s+Beta6(2);
M_l = 2/3*pi.*x(2).*(R1.^3-r^3);
M_s =  2/3*pi.*x(4).*(R2.^3-R1.^3);
M = mh+M_l+M_s;
v1 = (v0^2-(4.*ESE.*pi^3.*R2.^2.*x(3).^2.*x(3))./M).^0.5;
% Nonlinear constraints from first principles
c(1) = v1.*(M./(2.*pi.*R1.*sigma)).^0.5 - x(1);
c(2) = M_l+M_s-0.5;
c(3) = M_l.*EE_l+M_s.*EE_s-15e7; 
c(4) = M_l.*C_l+M_s.*C_s-10;
ceq = [];
end