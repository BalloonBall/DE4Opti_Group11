%% Objective function
% Minimising peak linear acceleration with variables:
% x(1)thickness of liner
% x(2)density of liner
% x(3)thickness of shell
% x(4)density of shell
function pla = PLA(x)
global Beta1 Beta4 v0 mh r
% Variable linking 
R1 = r+x(1);
R2 = r+x(1)+x(3);
sigma = Beta1(1).*x(2).^2+Beta1(2).*x(2)+Beta1(3);
ESE = Beta4(1).*x(4).^2+Beta4(2).*x(4)+Beta4(3);
M = mh+2/3*pi.*(x(2).*(R1.^3-r^3)+x(4).*(R2.^3-R1.^3));
v1 = (v0^2-(4.*ESE.*pi^3.*R2.^2.*x(3).^2.*x(3))./M).^0.5;
% Peak linear accleration from first principles
pla = v1./M.^0.5.*(2.*pi.*(R1).*(sigma)).^0.5/9.81;
end