function [x1,y1,aLo,dCLda] =  NACA_Airfoils(m,p,t,c,N,P)
x = linspace(0,c,(N+2)/2); % distance along chord x locations
for i = 1:length(x) % For each panel
    yt(i) = (t*c/0.2)*(0.2969*sqrt(x(i)/c) - 0.1260*(x(i)/c) - 0.3516*(x(i)/c)^2 + 0.2843*(x(i)/c)^3 - 0.1036*(x(i)/c)^4);
    if (0<x(i)) && (x(i)<=(p*c))
        yc(i) = m*(x(i)/p^2)*((2*p)-(x(i)/c));
    elseif (p*c<=x(i)) && (x(i)<=c)
        yc(i) = m*((c-x(i))/(1-p)^2)*(1+(x(i)/c)-(2*p));
    end
end
yc(1) = 0;
dx = x(2) - x(1);
dycdx = (gradient(yc))/dx;
xi = atan(dycdx);
xU = x - yt.*sin(xi);
xL = x + yt.*sin(xi);
yU = yc + yt.*cos(xi);
yL = yc - yt.*cos(xi);
x1 = [xU flip(xL(1:end-1))];
y1 = [yU flip(yL(1:end-1))];

if P ==1 % plot option
    figure
    plot(x1,y1)
    grid on
    axis('equal')
    title('NACA Airfoil Shape')
    ylabel('[m]')
    xlabel('[m]')
    hold off
end

ALPHA = linspace(-6,8,20);
VINF = 60;
for j = 1:length(ALPHA)
    [CL(j),~,~,~,~] = Vortex_Panel(x1,y1,VINF,ALPHA(j),0);
end
[h,~] = polyfit(-ALPHA,CL,1);
dCLda = h(1);
yint = h(2);
lobf = yint + -ALPHA.*dCLda;
aLo = -yint/dCLda;

if P == 1 %plot option
    figure
    hold on
    plot(-ALPHA,CL)
    plot(-ALPHA,lobf,'r+')
    plot(aLo,0,'b*')
    grid on
    title('NACA Airfoil Lift Slope')
    xlabel('Angle of Attack [alpha]')
    ylabel('Sectional Coefficient of Lift')
    legend('Angle of Attack vs Coefficient of Lift','Line of Best Fit','Zero Lift Angle of Attack','Location','SouthEast')
    hold off
end

end