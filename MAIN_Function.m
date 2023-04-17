% Alec Church
% ASEN 3111
% CA3 Flow over Finite Wings
% 4/3/22
% Code Assistance: Ben McHugh
clc
clear
close all
warning off
%% Problem 1
N = 100; %panels
c = 2; %chord length meters
P = 1; %plot bool

% NACA 0012
m = 0/100;
p = 0/10;
t = 12/100;
[x_0012,y_0012,aLo_0012,dCLda_0012] =  NACA_Airfoils(m,p,t,c,N,P);
fprintf('NACA 0012 Airfoil Zero Lift AoA: %f°.\n',aLo_0012);
fprintf('NACA 0012 Airfoil Lift Slope: %f \n',dCLda_0012);

% NACA 2412
m = 2/100;
p = 4/10;
t = 12/100;
[x_2412,y_2412,aLo_2412,dCLda_2412] =  NACA_Airfoils(m,p,t,c,N,P);
fprintf('NACA 2412 Airfoil Zero Lift AoA: %f°.\n',aLo_2412);
fprintf('NACA 2412 Airfoil Lift Slope: %f \n',dCLda_2412);

% NACA 4412
m = 4/100;
p = 4/10;
t = 12/100;
[x_4412,y_4412,aLo_4412,dCLda_4412] =  NACA_Airfoils(m,p,t,c,N,P);
fprintf('NACA 4412 Airfoil Zero Lift AoA: %f°.\n',aLo_4412);
fprintf('NACA 4412 Airfoil Lift Slope: %f \n',dCLda_4412);

%% Problem 2
clear  % managing workspace 
N = 50; %num odd terms for series expansion
P = 0;
b = 33 + 4/12; %span ft
c_t = 3 + 8.5/12; %tip chord
c_r = 5 + 4/12; %root chord
geo_t = 1;
geo_r = 0;
% NACA 0012
m = 0/100;
p = 0/10;
t = 12/100;
[x_tip,y_tip,alo_t,dclda_t] =  NACA_Airfoils(m,p,t,c_t,N,P);
[x_root,y_root,alo_r,dclda_r] =  NACA_Airfoils(m,p,t,c_r,N,P);

[e,c_L,c_Di] = PLLT(b,dclda_t,dclda_r,c_t,c_r,alo_t,alo_r,geo_t,geo_r,N);
V_inf = 82*1.687811;
rho_inf = 1.7556*10^-3;
L = (1/2)*c_L*rho_inf*V_inf^2*b*(c_t+c_r)/2;
D_i = (1/2)*c_Di*rho_inf*V_inf^2*b*(c_t+c_r)/2;
fprintf('The force of Lift is: %f\n',L)
fprintf('The force of Inducded Drag is: %f\n',D_i)

N10 = 2;
L0 = 0;
while abs(L0-L)>0.10
    [e,c_L,c_Di] = PLLT(b,dclda_t,dclda_r,c_t,c_r,alo_t,alo_r,geo_t,geo_r,N10);
    L0 = (1/2)*c_L*rho_inf*V_inf^2*b*(c_t+c_r)/2;
    N10 = N10+1;
end
fprintf('The number of odd terms required with 10 percent relative error: %d\n',N10)

N1 = 2;
L0 = 0;
while abs(L0-L)>0.01
    [e,c_L,c_Di] = PLLT(b,dclda_t,dclda_r,c_t,c_r,alo_t,alo_r,geo_t,geo_r,N1);
    L0 = (1/2)*c_L*rho_inf*V_inf^2*b*(c_t+c_r)/2;
    N1 = N1+1;
end
fprintf('The number of odd terms required with 1 percent relative error: %d\n',N1)

N01 = 2;
L0 = 0;
while abs(L0-L)>0.001
    [e,c_L,c_Di] = PLLT(b,dclda_t,dclda_r,c_t,c_r,alo_t,alo_r,geo_t,geo_r,N01);
    L0 = (1/2)*c_L*rho_inf*V_inf^2*b*(c_t+c_r)/2;
    N01 = N01+1;
end
fprintf('The number of odd terms required with 0.1 percent relative error: %d\n',N01)

%% Problem 3
j=0;
e_all = zeros(4,N);
for AR = [4,6,8,10]
    c_t1 = linspace(0,c_r,N)';
    c_r1 = c_r*ones(N,1);
    b = ((1/2)*AR.*(c_t1+c_r1))';
    j=j+1;
    for i = 1:N
       b_i = b(i);
       dclda_r = 2*pi;
       dclda_t = 2*pi;
       c_r_i = c_r1(i);
       c_t_i = c_t1(i);
       alo_r = 0.0011;
       alo_t = 0.0011;
       geo_r = 0.001;
       geo_t = 0.001;
       [e,c_L,c_Di] = PLLT(b_i,dclda_t,dclda_r,c_t_i,c_r_i,alo_t,alo_r,geo_t,geo_r,N01);
       e_all(j,i) = e;
    end
end
ctcr = c_t1./c_r1;
figure
hold on
plot(ctcr,e_all(1,:))
plot(ctcr,e_all(2,:))
plot(ctcr,e_all(3,:))
plot(ctcr,e_all(4,:))
title('Taper vs Span Efficiency Factor')
xlabel('Taper [Ct/Cr]')
ylabel('Span Efficiency Factor [e]')
legend('AR = 4','AR = 6','AR = 8','AR = 10','Location','SouthEast')