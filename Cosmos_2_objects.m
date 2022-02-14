clear; clc
% Jednostki:
%t - [doba]
%r - [1e6 km]
%m - [1e24 kg]

G = 4.9791e-04;     %Stała grawitacji w [(1e6 km)^3/(doba^2 * 1e24 kg)]
m1 = 2e6;           %Masa Słońca [1e24 kg]
m2 = 6;             %Masa Ziemi  [1e24 kg]
tmax = 1e3;         %Max. czas symulacji [dni]

x10 = 0;
y10 = 0;
x20 = 0;
y20 = 150;

v = 2.5732;      %Predkosc Ziemi na orbicie [1e6 km/doba]

vx10 = 0;        %Predkosc 'x' - Słońce
vy10 = 0;        %Predkosc 'y' - Słońce
vx20 = v;        %Predkosc 'x' - Ziemia
vy20 = 0;        %Predkosc 'y' - Ziemia

tspan = [0 tmax];
y0 = [x10; x20; y10; y20; vx10; vx20; vy10; vy20];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) rownania(t,y,m1,m2,G), tspan, y0, options);

n = 2e2;
Y = zeros(n,4);
tstep = 1/n:max(t)/n:max(t);
Y(:,1) = interp1(t,y(:,1),tstep);
Y(:,2) = interp1(t,y(:,2),tstep);
Y(:,3) = interp1(t,y(:,3),tstep);
Y(:,4) = interp1(t,y(:,4),tstep);


figure(1)
for i = 1: n
    plot(Y(i,1),Y(i,3),'r*',Y(i,2),Y(i,4),'bo');
    if (min(min(Y(:,1:2))) < max(max(Y(:,1:2))))...
    && (min(min(Y(:,3:4))) < max(max(Y(:,3:4))))
        xlim([min(min(Y(:,1:2))) max(max(Y(:,1:2)))]);
        ylim([min(min(Y(:,3:4))) max(max(Y(:,3:4)))]);
    end
    title({'Zmiana położenia w czasie obu obiektów';...
        strcat('Dzień nr:{ }',num2str(round(i.*max(t)./n)))})
    grid minor
    xlabel('x [mln km]')
    ylabel('y [mln km]')
% legend('Słońce','Ziemia')
%     axis equal
    movieVector(i) = getframe(figure(1));
end

%{
%Zapisanie filmu
myWriter = VideoWriter('Cosmos_2_2D.mp4');
myWriter.FrameRate = ceil(n/10);
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}
%{
figure(2)
plot(t,y(:,1),'.',t,y(:,2),'.');
title('Zmiana w czasie położeń ''x'' obiektów')
xlabel('t [days]')
ylabel('x [mln km]')
legend('Słońce','Ziemia')
%}

figure(3)
plot(y(:,1),y(:,3),'r*',y(:,2),y(:,4),'b.')
title('Rzut w czasie wszystkich położeń obu obiektów')
xlabel('x [mln km]')
ylabel('y [mln km]')
legend('Słońce','Ziemia')
axis equal



function [dydt] = rownania(t,y,m1,m2,G)
    r = sqrt((y(1)-y(2)).^2 + (y(3)-y(4)).^2);
    
    dydt = [y(5);...                 % dx1/dt
            y(6);...                 % dx2/dt
            y(7);...                 % dy1/dt
            y(8);...                 % dy2/dt
            (G.*m2./r.^2).*(y(2)-y(1))./r;...  % dvx1/dt
            (G.*m1./r.^2).*(y(1)-y(2))./r;...  % dvx2/dt
            (G.*m2./r.^2).*(y(4)-y(3))./r;...  % dvy1/dt
            (G.*m1./r.^2).*(y(3)-y(4))./r];    % dvy2/dt
end