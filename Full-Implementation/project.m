ma = 68; % kg
I = 0.38; % m
Ix = 12.1; % kg-m2
Iy = 11.67; % kg-m2
Iz = 1.08; % kg-m2
rw = 0.1; % m
Iw = 0.26; % kg-m2
mk = 19.6; % kg
rk = 0.22; % m
Ik = 0.38; % kg-m2
g = 9.98; % m/s2
alphaTemperature = (65.5/180)*pi; % degree temperature

Step1
% part a
syms yk(t) Xtheta(t) 
% Xtheta(t) = t definition
% yk(t) = t defnition
qx = transpose([yk(t) Xtheta(t)]);
qxDot(1,1) = diff(qx(1,1),t);
qxDot(2,1) = diff(qx(2,1),t);

syms t
Xr = transpose([transpose(qx) transpose(qxDot)]);

M11Qx = mk + Ik/(rk^2) + ma + (3*Iw*(cos(alphaTemperature)^2))/(2*rw*rw);
M12Qx = 3*Iw*(cos(alphaTemperature)^2)/(2*rw*rw)*rk - ma*I*cos(qx(2,1));
M21Qx = 3*Iw*(cos(alphaTemperature)^2)/(2*rw*rw)*rk - ma*I*cos(qx(2,1));
M22Qx =  ma*(I^2) + (3*Iw*(rk^2)*(cos(alphaTemperature)^2))/(2*(rw^2)) + Ix;
syms t
MQx(t) = [M11Qx M12Qx
        M21Qx M22Qx];
    
C11Qx = 0;
C12Qx = ma*I*qxDot(2,1)*sin(qx(2,1));
C21Qx = 0;
C22Qx = 0;
syms t
CQx(t) = [C11Qx C12Qx
      C21Qx C22Qx];

D11Qx = by*qxDot(1,1) + dy*sign(qxDot(1,1));
D12Qx = brx*qxDot(2,1) + drx*sign(qxDot(2,1));
syms t
DQx(t) = transpose([D11Qx D12Qx]);

G11Qx = 0;
G12Qx = (-1)*ma*g*I*sin(qx(2,1));

syms t
GQx(t) = transpose([G11Qx G12Qx]);
syms tauX(t)
% definition tau
Q11Qx = (1/rw)*tauX;
Q12Qx = (rk/rw)*tauX;
syms t
QQx(t) = transpose([Q11Qx Q12Qx]);
Xr(t) = transpose([transpose(qx) transpose(qxDot)]);
qxDotDot(t) = inv(MQx)*(QQx - CQx*qxDot - DQx - GQx); 
XrDot(t) = transpose([transpose(qxDot) transpose(qxDotDot)]);


% part b
syms xk(t) Ytheta(t) 
% Ytheta(t) = t definition
% xk(t) = t defnition
qy = transpose([xk(t) Ytheta(t)]);
qyDot(1,1) = diff(qy(1,1),t);
qyDot(2,1) = diff(qy(2,1),t);

syms t 
Xp = transpose([transpose(qy) transpose(qyDot)]);

M11Qy = mk + Ik/(rk^2) + ma + (3*Iw*(cos(alphaTemperature)^2))/(2*rw*rw);
M12Qy =  ma*I*cos(qy(2,1)) - (3*Iw*(cos(alphaTemperature)^2)/(2*rw*rw))*rk;
M21Qy = ma*I*cos(qy(2,1)) - (3*Iw*(cos(alphaTemperature)^2)/(2*rw*rw))*rk;
M22Qy =  ma*(I^2) + (3*Iw*(rk^2)*(cos(alphaTemperature)^2))/(2*(rw^2)) + Iy;

syms t
MQy(t) = [M11Qy M12Qy
          M21Qy M22Qy];
    
C11Qy = 0;
C12Qy = (-1)*ma*I*qxDot(2,1)*sin(qx(2,1));
C21Qy = 0;
C22Qy = 0;
syms t
CQy(t) = [C11Qy C12Qy
          C21Qy C22Qy];

D11Qy = bx*qyDot(1,1) + dx*sign(qyDot(1,1));
D12Qy = bry*qyDot(2,1) + dry*sign(qyDot(2,1));
syms t
DQy(t) = transpose([D11Qy D12Qy]);

G11Qy = 0;
G12Qy = (-1)*ma*g*I*sin(qy(2,1));
syms t
GQy(t) = transpose([G11Qy G12Qy]);
syms tauY(t)
% definition tau
Q11Qy = (-1)*(1/rw)*tauY;
Q12Qy = (rk/rw)*tauY;
syms t
QQy(t) = transpose([Q11Qy Q12Qy]);

Xp(t) = transpose([transpose(qy) transpose(qyDot)]);
qyDotDot(t) = inv(MQy)*(QQy - CQy*qyDot - DQy - GQy);
XpDot(t) = transpose([transpose(qyDot) transpose(qyDotDot)]);
%%% %%%
% part c
syms Ztheta(t) 
% Ztheta(t) = t definition
qz = Ztheta(t);
qzDot(1,1) = diff(qz,t);
syms t 
Xz(t) = transpose([qz qzDot]);
syms tauZ(t) t
qzDotDot(t) = (Ik/(Ik*Iz + 3*(Ik+Iz)*Iw*((rk^2)/(rw^2))*(sin(alphaTemperature)^2)))*(rk/rw)*tauZ;
syms t
XzDot(t) = transpose([qzDot qzDotDot]);

Step 2

% part A

% by = 0.01; % viscous damping
% brx = 0.01; % viscous damping
% dy = 0.01; % coulomb friction
% drx = 0.01; % coulomb friction

syms x y z m n by brx dy drx

% x = tauX
% y = yk
% z = thetaX
% m = ykDot
% n = thetaXDot

fr3_11 = (7.59*sin(z)-65.5*sin(2*z)+11.5*(n^2)*sin(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2) ;
fr3_21 = ((0.45)*by*m-0.03*brx*n+0.52*brx*n*cos(z))/((0.075*cos(z)-0.65*(cos(z)^2))+2.2);
fr3_31 = (0.45*dy*sign(m)-0.03*drx*sign(n))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2);
fr3_41 = (0.52*drx*sign(n)*cos(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2);
fr3(y,z,m,n,by,brx,dy,drx) = (-0.05)*(fr3_11 + fr3_21 + fr3_31 + fr3_41);
gr3(y,z,m,n,by,brx,dy,drx) = ((0.06*cos(z)+0.21)/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));
ykDotDotBlock(x,y,z,m,n,by,brx,dy,drx) = fr3 + gr3*x;
open_system('Step2A')
matlabFunctionBlock('Step2A/ykDotDotBlock',ykDotDotBlock)

syms x y z m n

% x = tauX
% y = yk
% z = thetaX
% m = ykDot
% n = thetaXDot

fr4_11 = (25.1*sin(z) + 0.037*(n^2)*sin(z) - 0.32*(z^2)*sin(2*z))/(0.075*cos(z) - 0.65*(cos(z)^2)+2.2);
fr4_21 = (-1)*((0.1*brx*n - 0.0015*by*m+ 0.025*by*m*cos(z))/(0.075*cos(z)-0.65*(cos(z^2))+2.2));
fr4_31 = (-1)*((0.1*drx*sign(n)-0.0015*dy*sign(m))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));
fr4_41 = ((0.025*dy*sign(m)*cos(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));

fr4(y,z,m,n,by,brx,dy,drx) = fr4_11 + fr4_21 + fr4_31 + fr4_41;
gr4(y,z,m,n,by,brx,dy,drx) = ((0.25*cos(z)+0.21)/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));
thetaXDotDotBlock(x,y,z,m,n,by,brx,dy,drx) = fr4 + gr4*x;

open_system('Step2A')
matlabFunctionBlock('Step2A/thetaXDotDotBlock',thetaXDotDotBlock)

% part B

% bx = 0 % 0.01; % viscous damping
% bry = 0 %0.01; % viscous damping
% dx = 0 %0.01; % coulomb friction
% dry = 0 %0.01; % coulomb friction

syms x y z m n bx bry dx dry

% x = tauY
% y = xk
% z = thetaY
% m = xkDot
% n = thetaYDot

fp3_11 = (65.5*sin(2*z)-7.59*sin(z)-11.3*(n^2)*sin(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2) ;
fp3_21 = (0.03*bry*n + 0.44*bx*m-0.52*bry*n*cos(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2);
fp3_31 = (0.03*dry*sign(n)+0.44*dx*sign(m))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2);
fp3_41 = (-1)*(0.52*dry*sign(n)*cos(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2);
fp3(y,z,m,n,bx,bry,dx,dry) = (-0.05)*(fp3_11 + fp3_21 + fp3_31 - fp3_41);
gp3(y,z,m,n,bx,bry,dx,dry) = (-1)*((0.06*cos(z)+0.21)/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));
xkDotDotBlock(x,y,z,m,n,bx,bry,dx,dry) = fp3 + gp3*x;
open_system('Step2B')
matlabFunctionBlock('Step2B/xkDotDotBlock',xkDotDotBlock)

syms x y z m n

% x = tauY
% y = xk
% z = thetaY
% m = xkDot
% n = thetaYDot

fp4_11 = (-25.1*sin(z) - 0.037*(n^2)*sin(z) + 0.32*(z^2)*sin(2*z))/(0.075*cos(z) - 0.65*(cos(z)^2)+2.2);
fp4_21 = ((0.1*bry*n + 0.0015*bx*m - 0.025*bx*m*cos(z))/(0.075*cos(z)-0.65*(cos(z^2))+2.2));
fp4_31 = ((0.1*dry*sign(n) + 0.0015*dx*sign(m))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));
fp4_41 = (-1)*((0.025*dx*sign(m)*cos(z))/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));

fp4(y,z,m,n,bx,bry,dx,dry) = (-1)*(fp4_11 + fp4_21 + fp4_31 + fp4_41);
gp4(y,z,m,n,bx,bry,dx,dry) = ((0.25*cos(z)+0.21)/(0.075*cos(z)-0.65*(cos(z)^2)+2.2));
thetaYDotDotBlock(x,y,z,m,n,bx,bry,dx,dry) = fp4 + gp4*x;

open_system('Step2B')
matlabFunctionBlock('Step2B/thetaYDotDotBlock',thetaYDotDotBlock)
  

Step 3
syms alpha
torque1(alpha) = [(2/(3*cos(alpha))) 0 1/(3*sin(alpha))];
torque2(alpha) = [(-1/(3*cos(alpha))) sqrt(3)/(3*cos(alpha)) 1/(3*sin(alpha))];
torque3(alpha) = [(-1/(3*cos(alpha))) ((-1)*sqrt(3))/(3*cos(alpha)) 1/(3*sin(alpha))];
torqueConversion = [torque1;torque2;torque3];
syms tauX tauY tauZ Tau
Tau = [tauX;tauY;tauZ];

conversedTau(alpha,Tau(1,1),Tau(2,1),Tau(3,1)) = torqueConversion*Tau;

% considering alpha = 66.5
consideredAlpha = (66.5/180)*pi;
tau = conversedTau(consideredAlpha,Tau(1,1),Tau(2,1),Tau(3,1));
tau1 = tau(1,1);
tau2 = tau(2,1);
tau3 = tau(3,1);

open_system('Step3');
matlabFunctionBlock('Step3/torqueConversion',tau1,tau2,tau3);

Step  4
syms thetaZDot theta thetaZDestination kdz kpz
TauZ(thetaZDot,theta) = kdz*thetaZDot + kpz*(theta - thetaZDestination);
open_system('Step4A');
matlabFunctionBlock('Step4A/PIDController',TauZ);

Step 5
LQR controller
Kr = [-6.98 -4.45 20.88 6.09 5.11];
Kp = [-3.27  4.01 29.07 8.54 5.72];
syms thetaXDot thetaX ykDot yk ykd x5r
vy(thetaXDot,thetaX,ykDot,yk,ykd,x5r) = (-1)*Kr*transpose([thetaXDot thetaX ykDot yk-ykd x5r]);
syms thetaYDot thetaY xkDot xk xkd x5p
vx(thetaYDot,thetaY,xkDot,xk,xkd,x5p) = (-1)*Kp*transpose([thetaYDot thetaY xkDot xk-xkd x5p]);
open_system('LQRController.slx');
matlabFunctionBlock('LQRController/LQR',vx,vy);

PI controller
syms kessi Wn by brx dy drx
syms z % considered z as thetaX
gr3inverse(z,by,brx,dy,drx) = 1/gr3(0,z,0,0,by,brx,dy,drx);
syms bx bry dx dry
gp3inverse(z,bx,bry,dx,dry) = 1/gp3(0,z,0,0,bx,bry,dx,dry);

% doesnt match
%figure
%fplot(gp3inverse);
%xlim([1 10]);
%ylim([0 10]);
%title('gp3Inverse');

zero = 0;
Kpr(kessi,Wn,by,brx,dy,drx) = 2*gr3inverse(zero,by,brx,dy,drx)*kessi*Wn;
Kir(kessi,Wn,by,brx,dy,drx) = gr3inverse(zero,by,brx,dy,drx)*(Wn^2);

zero = 0.001;
Kpp(kessi,Wn,bx,bry,dx,dry) = 2*gp3inverse(zero,bx,bry,dx,dry)*kessi*Wn;
Kip(kessi,Wn,bx,bry,dx,dry) = gp3inverse(zero,bx,bry,dx,dry)*(Wn^2);

syms vx vy xk yk xkDot ykDot x5p x5r by brx dy drx bx bry dx dry
uPIy(vx,xkDot,x5p,xk,bx,bry,dx,dry) = Kpp*(vx - xkDot) + Kip*(x5p - xk);
uPIx(vy,ykDot,x5r,yk,by,brx,dy,drx) = Kpr*(vy - ykDot) + Kir*(x5r - yk);

open_system('PIController.slx');
matlabFunctionBlock('PIController/PIController_KP_KI',uPIx,uPIy);

Feedforward Compensation Term
zero = 0.000001;
syms by brx dy drx
Uxf = vpa((-1)*gr3inverse(0,by,brx,dy,drx)*fr3(zero,zero,zero,zero,by,brx,dy,drx));
syms bx bry dx dry
Uyf = vpa((-1)*gp3inverse(0,bx,bry,dx,dry)*fp3(zero,zero,zero,zero,bx,bry,dx,dry));
open_system('FeedforWardCompensationTerm.slx');
matlabFunctionBlock('FeedforWardCompensationTerm/frictions',Uxf,Uyf);



