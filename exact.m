%%    Report Joukowski airfoil
%     
%  Name:- Mohamed Abd El Mawgoud Ghomeam
%  Sec:- 2     B.N:- 16
%     
%      
%% useful command
close all
clearvars
clc

%% inputs
cord=1;                  % Chord length of the airfoil (m)
Vinf=100;                  % Free stream velocity ( m/s)
maxthik=.07;               % max. thickness to chord ratio
maxcamb=.05;               % max. camber to chord ratio 
alfad = 8;                 % Angle of attack  (deg)

%% important calculations

b=cord/4;
e=maxthik/1.3;
beta=2*maxcamb;
a=b*(1+e)/cos(beta);
alfa = alfad * pi / 180; 
xo = -b * e;
yo = a * sin(beta);
theta_=0:pi/100:2*pi;
r_=a;
theta=atan2(r_*sin(theta_)+a*beta,r_*cos(theta_)-b*e);
r=b*(1+e*(1-cos(theta))+beta*sin(theta));



%% airfoil coordinate x1&y1
theta_=0:pi/100:2*pi;
x1_app=2*b*cos(theta_);
y1_app=2*b*e*(1-cos(theta_)).*sin(theta_)+2*b*beta*sin(theta_).^2;
sizeuplo=size(x1_app,2)-1;
x1_up=x1_app(1:sizeuplo/2+1); %X-coord for upper surface Airfoil
y1_up=y1_app(1:sizeuplo/2+1); %Y-coord for upper surface Airfoil
x1_lo=x1_app(sizeuplo:-1:sizeuplo/2); %X-coord for lower surface Airfoil
y1_lo=y1_app(sizeuplo:-1:sizeuplo/2); %Y-coord for lower surface Airfoil
y1_camb=(y1_lo+y1_up)/2; %y-coord for camber line
figure
plot(x1_up,y1_camb,'r--','LineWidth',1) ; hold on %color red camber line
plot(x1_up,y1_up,'g-','LineWidth',1) ; hold on  %color green for upper
plot(x1_lo,y1_lo,'b-','LineWidth',1)  %color blue for lower
xlabel('$X_{B}$ $for$ $Airfoil$','interpreter','latex','FontSize',14);
ylabel('$Y_{B}$ $for$ $Airfoil$','interpreter','latex','FontSize',14);
title('$Airfoil$ $shape$ $approximate$','interpreter','latex','FontSize',14);
axis equal
grid on

%% streamlines (note the airfoil is at angle of attac = alfad)
theta__=0:pi/100:2*pi;
r__=a:a/80:6*a; 
[theta__,r__]=meshgrid(theta__,r__);
alpha__=alfad*pi/180;
x_=r__.*cos(theta__);
y_=r__.*sin(theta__);
x=x_+xo;
y=y_+yo;
x1=x.*(1+b^2./(x.^2+y.^2));
y1=y.*(1-b^2./(x.^2+y.^2));
epsi=Vinf*(r__.*sin(theta__-alpha__)+a^2./r__.*sin(alpha__-theta__)...
    +2*a*sin(alpha__+beta).*log(r__));
%transform
%Here I rotate airfoil and keep flow horizontal
xx1=x1*cos(alpha__)+y1*sin(alpha__);
yy1=-x1*sin(alpha__)+y1*cos(alpha__);
figure
fill(xx1(1,1:end),yy1(1,1:end),'r');   hold on;     grid on
contour(xx1,yy1,epsi,76,'b','LineWidth',1);
axis([-1 1 -0.6 0.6]);
xlabel('$airfoil$ $X-axis$','interpreter','latex','FontSize',14)
ylabel('$airfoil$ $Y-axis$','interpreter','latex','FontSize',14)
title('$Streamlines$ $around$ $the$ $airfoil.$','interpreter',...
    'latex','FontSize',14);
grid on

%% V/Vinf versus x1 & cp versus xi
clear vt_ V
vt_=-Vinf*(2*sin(theta_-alfa)+2*sin(alfa+beta));
V=sqrt(vt_.^2./(1-2*(b^2./r.^2).*cos(2*theta)+(b^4./r.^4)));
cp=1-V.^2/Vinf^2;
sizeuplo=size(x1_app,2)-1;
V_up=V(1:sizeuplo/2+1);
cp_up=cp(1:sizeuplo/2+1);
V_lo=V(sizeuplo:-1:sizeuplo/2);
cp_lo=cp(sizeuplo:-1:sizeuplo/2);

figure
plot(x1_up,cp_up,'g-','LineWidth',1)  ;  hold on %color green for upper
plot(x1_lo,cp_lo,'b-','LineWidth',1)   %color blue for lower
xlabel('$X_{B}$','interpreter','latex','FontSize',14);
ylabel('$C_{p}$','interpreter','latex','FontSize',14);
title('$Pressure$ $over$ $airfoil$','interpreter','latex','FontSize',14);
grid on

figure
plot(x1_up,V_up/Vinf,'g-','LineWidth',1) ; hold on %color green for upper
plot(x1_lo,V_lo/Vinf,'b-','LineWidth',1)  %color blue for lower
xlabel('$X_{B}$','interpreter','latex','FontSize',14);
ylabel('$\frac{V}{V_{inf}}$','interpreter','latex','FontSize',14);
title('$Velocity$ $over$ $airfoil$','interpreter','latex','FontSize',14);
grid on

%% lift coeff Cl Vs alpha
figure
alpha_cl=(-10:0.1:10)*pi/180;
Cl=2*pi.*(1+e).*sin(alpha_cl+beta);
plot(alpha_cl/pi*180,Cl,'g-','LineWidth',2)
xlabel('$alpha$','interpreter','latex','FontSize',14);
ylabel('$C_{L}$','interpreter','latex','FontSize',14);
title('$C_{L}$ $approximate$','interpreter','latex','FontSize',14);
grid on

%% Cm calculation Vs alpha
alpha_=(-10:.1:10)*pi/180;
Cp=zeros(1,size(theta_,2));
Cm=zeros(1,size(alpha_,2));

for i=1:size(alpha_,2)
    vt_=-Vinf*(2*sin(theta_-alpha_(i))+2*sin(alpha_(i)+beta));
    V=sqrt(vt_.^2./(1-2*(b^2./r.^2).*cos(2*theta)+(b^4./r.^4)));
    Cp=(1-V.^2/Vinf.^2);
    M=0;
    %I will show some rules in outer paper help me to get the relation
    %I get Cm with numerical method
   for j=1:size(theta_,2)-1 
    M=M-(Cp(j)+Cp(j+1))/2*...
        ((x1_app(j)+x1_app(j+1))/2*(x1_app(j+1)-x1_app(j)));
   end
   Cm(i)=M/(cord^2);
end

figure;
plot(alpha_*180/pi,Cm,'r-','LineWidth',2);
xlabel('$alpha$','interpreter','latex','FontSize',14);
ylabel('$C_{m}$','interpreter','latex','FontSize',14);
title('$C_{m}$ $over$ $airfoil$ $using$ $numerical$ $method$',...
    'interpreter','latex','FontSize',14);
grid on;