clear; clc
% Jednostki:
%t - [doba]
%r - [1e6 km]
%m - [1e24 kg]

%% Dane (do zmieniania wedle woli)
G = 4.9791e-04;   %Stała grawitacji w [(1e6 km)^3/(doba^2 * 1e24 kg)]
% W rzeczywistości: G = 4.9791e-04 [jednostki j.w.]
m1 = 2e4;               %Masa Ciała 1 (np. Słońca) [1e24 kg]
m2 = 2e4;               %Masa Ciała 2 (np. Ziemi)  [1e24 kg]
m3 = 2e6;               %Masa Ciała 3              [1e24 kg]
tmax = 6.6e5;             %Max. czas symulacji [dni]

x10 = 0;                %[mln km] - Początkowe położenie ciała 1 na 'x'
y10 = 0;                %[mln km] - Początkowe położenie ciała 1 na 'y'
x20 = 1e2;              %[mln km] - Początkowe położenie ciała 2 na 'x'
y20 = 0;                %[mln km] - Początkowe położenie ciała 2 na 'y'
x30 = 5e2;              %[mln km] - Początkowe położenie ciała 3 na 'x'
y30 = 5e2;              %[mln km] - Początkowe położenie ciała 3 na 'y'

v = 2.5732e-2;          %Predkosc Ziemi na orbicie [mln km/doba]

vx20 = 0;               %[mln km/dzien] Predkosc początkowa 'x' - Ciało 2 np. Ziemia
vy20 = 0;               %[mln km/dzien] Predkosc początkowa 'y' - Ciało 2 np. Ziemia
vx30 = 0;               %[mln km/dzien] Predkosc początkowa 'x' - Ciało 3
vy30 = 0;               %[mln km/dzien] Predkosc początkowa 'y' - Ciało 3
vx10 = 0;
vy10 = 0;

% Prawo zachowania pędu, by się gwiazda nie ruszała
% vx10 = (-m2.*vx20-m3.*vx30)./m1;  %[mln km/dzien] Predkosc początkowa 'x' - Ciało 1 np. Słońce
% vy10 = (-m2.*vy20-m3.*vy30)./m1;  %[mln km/dzien] Predkosc początkowa 'y' - Ciało 1 np. Słońce

%% Obliczenia
tspan = [0 tmax];
y0 = [x10;  y10;  x20;  y20;  x30;  y30;...
      vx10; vy10; vx20; vy20; vx30; vy30];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) rownania(t,y,m1,m2,m3,G), tspan, y0, options);

n = 3e2;
Y = zeros(n,6);
tstep = 1/n:max(t)/n:max(t);
Y(:,1) = interp1(t,y(:,1),tstep);
Y(:,2) = interp1(t,y(:,2),tstep);
Y(:,3) = interp1(t,y(:,3),tstep);
Y(:,4) = interp1(t,y(:,4),tstep);
Y(:,5) = interp1(t,y(:,5),tstep);
Y(:,6) = interp1(t,y(:,6),tstep);

%Wprowadzić ślad ruchu do wykresu!
figure(1)
set(0, 'defaultfigureposition', [500, 0, 1000, 1000]);
for i = 1: n
    plot(Y(i,1),Y(i,2),'r*',Y(i,3),Y(i,4),'co',Y(i,5),Y(i,6),'go');
    xlim([min(min([Y(:,1),Y(:,3),Y(:,5)])) max(max(([Y(:,1),Y(:,3),Y(:,5)])))]);
    ylim([min(min(([Y(:,2),Y(:,4),Y(:,6)]))) max(max(([Y(:,2),Y(:,4),Y(:,6)])))]);
    title({'Zmiana położenia w czasie obiektów';...
        strcat('Dzień nr:{ }',num2str(round(i.*max(t)./n)))})
    grid on
	set(gca,'Color',[0.1 0.1 0.1]);     %background
    ax = gca;
    ax.GridColor = [0.8, 0.8, 0.8];     %gridlines
    xlabel('x [mln km]')
    ylabel('y [mln km]')
% legend('Ciało 1','Ciało 2','Ciało 3')
%     axis equal
    movieVector(i) = getframe(figure(1));
end


%Zapisanie filmu
myWriter = VideoWriter('Cosmos_3_2D.mp4');
myWriter.FrameRate = ceil(n/10);
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}

figure(2)
plot(t,y(:,1),'r-',t,y(:,3),'b-',t,y(:,5),'g-');
title('Zmiana w czasie położeń ''x'' obiektów')
xlabel('t [days]')
ylabel('x [mln km]')
legend('Ciało 1','Ciało 2','Ciało 3')
%}

figure(3)
plot(y(:,1),y(:,2),'r.',y(:,3),y(:,4),'b-',y(:,5),y(:,6),'g-')
title('Tory ruchu wszystkich obiektów')
xlabel('x [mln km]')
ylabel('y [mln km]')
legend('Ciało 1','Ciało 2','Ciało 3')
axis equal

%% Funkcja z równaniem różniczkowym

% Wprowadzić teorię względności, bo przekroczyło mi prędkość światła!
% Albo napisać osobną funkcję na to...

function [dydt] = rownania(t,y,m1,m2,m3,G)
    r12 = sqrt((y(1)-y(3)).^2 + (y(2)-y(4)).^2);
    r13 = sqrt((y(1)-y(5)).^2 + (y(2)-y(6)).^2);
    r23 = sqrt((y(3)-y(5)).^2 + (y(4)-y(6)).^2);
    
    dydt = [y(7);...                 % dx1/dt
            y(8);...                 % dy1/dt
            y(9);...                 % dx2/dt
            y(10);...                % dy2/dt
            y(11);...                % dx3/dt
            y(12);...                % dy3/dt
            
            G.*(m2.*(y(3)-y(1))./(r12.^3) + m3.*(y(5)-y(1))./(r13.^3));...  % dv1x/dt
            G.*(m2.*(y(4)-y(2))./(r12.^3) + m3.*(y(6)-y(2))./(r13.^3));...  % dv1y/dt
            
            G.*(m1.*(y(1)-y(3))./(r12.^3) + m3.*(y(5)-y(3))./(r23.^3));...  % dv2x/dt
            G.*(m1.*(y(2)-y(4))./(r12.^3) + m3.*(y(6)-y(4))./(r23.^3));...  % dv2y/dt
            
            G.*(m1.*(y(1)-y(5))./(r13.^3) + m2.*(y(3)-y(5))./(r23.^3));...  % dv3x/dt
            G.*(m1.*(y(2)-y(6))./(r13.^3) + m2.*(y(4)-y(6))./(r23.^3))];    % dv3y/dt
        
end



