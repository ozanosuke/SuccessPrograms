function step
clear all
global s1 s2 l1 l2 m1 m2 I1 I2 g v1 cx1 cy1 v2 cx2 cy2

l1=1;
l2=1;
s1=0.5;
s2=0.5;
g=9.81;
m1=1;
m2=1;
I1=m1*l1^2/12;
I2=m2*l2^2/12;

theta1=pi*2/3;
theta2=pi/3;

TSPAN=[0:1E-2:16];
X0=[s1*cos(theta1) s1*sin(theta1) theta1 l1*cos(theta1)-(l1-s1)*cos(theta2) l1*sin(theta1)-(l1-s1)*sin(theta2) theta2 0 0 0 0 0 0]';
[t,X]=ode45(@func,TSPAN,X0);

for i=1:length(X)
    figure(1);
    plot([0,0],[0,0]);
    grid;
    axis square;
    xmin=-2; xmax=4; ymin=-2; ymax=4;
    axis([xmin,xmax,ymin,ymax]);
    hold on;
    xQ=X(i,1)-s1*cos(X(i,3));
    yQ=X(i,2)-s1*sin(X(i,3));
    xP1=X(i,1)+(l1-s1)*cos(X(i,3));
    yP1=X(i,2)+(l1-s1)*sin(X(i,3));
    xP2=X(i,4)+(l2-s2)*cos(X(i,6));
    yP2=X(i,5)+(l2-s2)*sin(X(i,6));
    xR=X(i,4)-s2*cos(X(i,6));
    yR=X(i,5)-s2*sin(X(i,6));
    plot([xQ,xP1],[yQ,yP1],'linewidth',4);
    plot([xR,xP2],[yR,yP2],'linewidth',4);
    drawnow
    hold off
    set(1,'doublebuffer','on')
end

function dXdt=func(t,X)
global s1 s2 l1 l2 m1 m2 I1 I2 g v1 cx1 cy1 v2 cx2 cy2 theta1 theta2

theta1=pi*2/3;
theta2=pi/3;
if 0<=t&&t<=4
    v1=0;
    v2=0.5;
    cx1=0;
    cx2=-1;
    cy1=0;
    cy2=0;
elseif 4<t&&t<=8
    v1=0.5;
    v2=0;
    cx1=-2;
    cx2=1;
    cy1=0;
    cy2=0;
elseif 8<t&&t<=12
    v1=0;
    v2=0.5;
    cx1=2;
    cx2=-3;
    cy1=0;
    cy2=0;
elseif 12<t&&t<=16
    v1=0.5;
    v2=0;
    cx1=-4;
    cx2=3;
    cy1=0;
    cy2=0;    
end


M1=[m1 0 0;
    0 m1 0;
    0 0 I1];
M2=[m2 0 0;
    0 m2 0;
    0 0 I2];
Q=[0;
    -m1*g;
    0;
    0;
    -m2*g;
    0];
Cq1=[1 0 -(l1-s1)*sin(X(3));
    0 1 (l1-s1)*cos(X(3));
    1 0 s1*sin(X(3));
    0 1 -s1*cos(X(3));
    0 0 0;
    0 0 0];

Cq2=[-1 0 (l2-s2)*sin(X(6));
    0 -1 -(l2-s2)*cos(X(6));
    0 0 0;
    0 0 0;
    1 0 s2*sin(X(6));
    0 1 -s2*cos(X(6))];
A=[M1         zeros(3,3) Cq1';
   zeros(3,3) M2         Cq2';
   Cq1        Cq2        zeros(6,6)];
 
alph=500;
beta=500;

C=[X(1)+(l1-s1)*cos(X(3))-X(4)-(l2-s2)*cos(X(6));
    X(2)+(l1-s1)*sin(X(3))-X(5)-(l2-s2)*sin(X(6));
    X(1)-s1*cos(X(3))-(v1*t+cx1);
    X(2)-s1*sin(X(3))-cy1;
    X(4)-s2*cos(X(6))-(v2*t+cx2);
    X(5)-s2*sin(X(6))-cy2];

C1=[X(7)-(l1-s1)*sin(X(3))*X(9)-X(10)+(l2-s2)*sin(X(6))*X(12);
    X(8)+(l1-s1)*cos(X(3))*X(9)-X(11)-(l2-s2)*cos(X(6))*X(12);
    X(7)+s1*sin(X(3))*X(9)-v1;
    X(8)-s1*cos(X(3))*X(9);
    X(10)+s2*sin(X(6))*X(12)-v2;
    X(11)-s2*cos(X(6))*X(12)];

Gm=[-(l1-s1)*cos(X(3))*X(9)^2+(l2-s2)*cos(X(6))*X(12)^2;
    -(l1-s1)*sin(X(3))*X(9)^2+(l2-s2)*sin(X(6))*X(12)^2;
    s1*cos(X(3))*X(9)^2;
    s1*sin(X(3))*X(9)^2;
    s2*cos(X(6))*X(12)^2;
    s2*sin(X(6))*X(12)^2]-2*alph*C1-(beta^2)*C;

RHS=[Q;
    Gm];
ACC=A\RHS;
dXdt=zeros(6,1);
dXdt(1:6)=X(7:12);
dXdt(7:12)=ACC(1:6);
