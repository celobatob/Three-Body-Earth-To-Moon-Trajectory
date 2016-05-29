%*** 3-body Earth-Moon Problem *****************************************
clear all

%*** Earth Frame *******************************************************
Eradius = 6371; %km
muE = 3.986e5;

%*** LEO ********************************
e_leo = 0;
alt_leo = 322; %km

r_leo = Eradius + alt_leo; %km
v_leo = sqrt(muE/r_leo);

%v0 = v_leo; % For LEO

% Moon Trajectory Initial Conditions
%********************
 v0 = 10.78690; %km/s
 gma0 = 132.07500; %deg
% v0 = 10.78590;
%gma0 = 142.47500;
%********************

%*** Frame Conversion ****************************************************
rM = 384400; %km
wM = 2.6491e-6; %rad/s
mu = 0.01215;

DU = rM;
VU = wM*DU;
TU = DU/VU;
WU = 1;

%*** Convert to Rotating Frame ************************
% Position Components (Rotating Frame)
x0 = (r_leo*cosd(180-gma0))/rM + mu; %DU
y0 = (r_leo*sind(180-gma0))/rM; %DU

% Position Rotating Frame
r0s = [x0 y0 0]

% Velocity Components (Earth Frame)
vx0 = (v0*cosd(270-gma0))/VU; %DU/TU
vy0 = (v0*sind(270-gma0))/VU; %DU/TU

% Velocity Rotating Frame
v0s = [vx0 vy0 0] - cross([0 0 1],r0s)

%*** ODE Solver **************************************
Ini = [x0 v0s(1) y0 v0s(2)]
tf = 5;
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Refine',10);
[T,Y] = ode45(@func,[0 tf],Ini,options);

%*** Plot System **************************************
% Earth
deg = 0:360;
xE = (Eradius/rM)*cosd(deg) + mu;
yE = (Eradius/rM)*sind(deg);
figure(1)
plot(xE,yE,'c','LineWidth',2)
hold on

% Moon
Mradius = 1738; %km
xM = (Mradius/rM)*cosd(deg) - (1-mu);
yM = (Mradius/rM)*sind(deg);
plot(xM,yM,'m','LineWidth',2)

% Orbit Trajectory
plot(Y(:,1),Y(:,3),'b')
title('Rotating Frame')
xlabel('x (DU)')
ylabel('y (DU)')
axis equal
hold off

%*** Mission Objective Check *********************************************
% Periselnium Altitude
Dm = sqrt((Y(:,1) + (1-mu)).^2 + Y(:,3).^2); %DU
[r_per,ind] = min(Dm); %DU
per_rad = r_per*rM
per_alt = per_rad - Mradius %km
disp('km')

% Earth Altitude
n = numel(T);
De = sqrt((Y(10000:n,1) - mu).^2 + Y(10000:n,3).^2); %DU
rE_per = min(De); %DU
E_alt = rE_per*rM - Eradius %km
disp('km')

%*** Plot in Earth Inertial Frame ****************************************
% Convert to IJK
rS = zeros(3,n);
rS(1,:) = Y(:,1)-mu;
rS(2,:) = Y(:,3);
for i = 1:n
    Gam = WU*T(i);
    C = [cos(Gam) sin(Gam) 0; -sin(Gam) cos(Gam) 0; 0 0 1];
    rijk(:,i) = inv(C)*rS(:,i);
end

% Earth
xE = (Eradius/rM)*cosd(deg);
yE = (Eradius/rM)*sind(deg);
figure(2)
plot(xE,yE,'c','LineWidth',2)
hold on

% S/C Trajectory

plot(rijk(1,:),rijk(2,:))
title('Earth Inertial Frame')
xlabel('x (DU)')
ylabel('y (DU)')

% Moon Trajectory
theta = WU*T;
xMoon = -cos(theta);
yMoon = -sin(theta);
plot(xMoon,yMoon,'r')


% Moon
thetaM = WU*T(8688);
xmMoon = -cos(thetaM);
ymMoon = -sin(thetaM);
xM = (Mradius/rM)*cosd(deg) + xmMoon;
yM = (Mradius/rM)*sind(deg) + ymMoon;
plot(xM,yM,'m','LineWidth',2)
axis equal
hold off


