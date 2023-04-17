function [e,c_L,c_Di] = PLLT(b,dclda_t,dclda_r,c_t,c_r,alo_t,alo_r,geo_t,geo_r,N)
d = 1:N;
theta = (pi/(2*N))*d';
S = b*(c_t+c_r)/2;
AR = b^2/S;
a = (dclda_t-dclda_r)*cos(theta)+dclda_r;
c = (c_t-c_r)*cos(theta)+c_r;
aGeo = (geo_t-geo_r)*cos(theta)+geo_r;
a0 = (alo_t-alo_r)*cos(theta)+alo_r;

for i = 1:N
    alpha_eff(i) = aGeo(i) - a0(i);
    for j = 1:N
        f(i,j) = 4*b/(a(i)*c(i))*sin((2*j-1)*theta(i))+(2*j-1)*sin((2*j-1)*theta(i))/sin(theta(i));
    end
end

A = linsolve(f,alpha_eff');
d = 0;
for k = 2:N
        d=d+(2*k-1)*(A(k)/A(1))^2;
end
e = (1+d)^-1;
c_L = A(1)*pi*AR;
c_Di = c_L^2/(pi*e*AR);
end