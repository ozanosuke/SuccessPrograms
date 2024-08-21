function Example5_3
clear all
global m1 m2 I1 I2 g l1 l2 s1 s2 v cx cy
m1=1;
m2=1;
l1=1;
l2=1;
s1=0.5;
s2=0.5;
g=9.81;
I1=m1*l1^2/12;
I2=m2*l2^2/12;
v=0.5;
cx=l1*cos(pi/4)+l2*cos(pi*2/3);
cy=l1*sin(pi/4)+l2*sin(pi*2/3);

TSPAN=[0:1E-2:2];
X0=[s1*cos(pi/4) s1*sin(pi/4) pi/4 l1*cos(pi/4)+s2*cos(pi*2/3) l1*sin(pi/4)+s2*sin(pi*2/3) pi*2/3 0 0 0 0 0 0];
[t,X]=ode45(@func,TSPAN,X0);

for i=1:length(X)
    figure(1);
    plot([0,0],[0,0]);
    grid;
    axis square;
    xmin=-1;xmax=2;ymin=0;ymax=3;
    axis([xmin,xmax,ymin,ymax]);
    hold on;
    xQ=X(i,1)-s1*cos(X(i,3));
    yQ=X(i,2)-s1*sin(X(i,3));
    xP1=X(i,1)+(l1-s1)*cos(X(i,3));
    yP1=X(i,2)+(l1-s1)*sin(X(i,3));
    xP2=X(i,4)-s2*cos(X(i,6));
    yP2=X(i,5)-s2*sin(X(i,6));
    xR=X(i,4)+(l2-s2)*cos(X(i,6));
    yR=X(i,5)+(l2-s2)*sin(X(i,6));
    
    plot([xQ,xP1],[yQ,yP1],'linewidth',4);
    plot([xR,xP2],[yR,yP2],'linewidth',4);
    drawnow
    hold off
    set(1,'doublebuffer','on')
end

figure(2);
plot(t,X(:,5)+(l2-s2)*sin(X(:,6)))
xlabel('time','FontSize',18);
ylabel('Ry','FontSize',18);

function dXdt=func(t,X)
global m1 m2 I1 I2 g l1 l2 s1 s2 v cx cy
M=[m1 0 0 0 0 0;
    0 m1 0 0 0 0;
    0 0 I1 0 0 0;
    0 0 0 m2 0 0;
    0 0 0 0 m2 0;
    0 0 0 0 0 I2];
Qg=[0;
    -m1*g;
    0;
    0;
    -m2*g;
    0];
Qf=[0;
    0;
    (I1*((l2*v^2*cos(X(3))^2*cos(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (v^2*cos(X(3))*cos(X(6))^3)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))) - (I2*((v^2*cos(X(3))^2*cos(X(6))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^3*cos(X(6)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(3))*sin(X(6)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + g*l1*m2*cos(X(3)) + g*m1*s1*cos(X(3)) + g*m2*s2*cos(X(6)) + (l1*m2*cos(X(3))*(l2 - s2)*((v^2*cos(X(3))^2*cos(X(6))^3)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^3*sin(X(3))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^3*cos(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^3*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + (m2*s2*cos(X(6))*(l2 - s2)*((v^2*cos(X(3))^2*cos(X(6))^3)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^3*sin(X(3))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^3*cos(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^3*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + (m1*s1^2*sin(X(3))*((v^2*cos(X(3))^2*cos(X(6))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^2*cos(X(6))^2*sin(X(3)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(3))*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))) + (l1*m2*sin(X(3))*(l2 - s2)*((v^2*cos(X(3))^2*cos(X(6))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^2*cos(X(6))^2*sin(X(3)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(3))*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + (m2*s2*sin(X(6))*(l2 - s2)*((v^2*cos(X(3))^2*cos(X(6))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^2*cos(X(6))^2*sin(X(3)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(3))*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + (m1*s1^2*cos(X(3))*((v^2*cos(X(3))^2*cos(X(6))^3)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^3*sin(X(3))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^3*cos(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^3*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3))));
    0;
    0;
    g*m2*s2*cos(X(6)) - (I2*((v^2*cos(X(3))^2*cos(X(6))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^3*cos(X(6)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(3))*sin(X(6)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + (m2*s2*cos(X(6))*(l2 - s2)*((v^2*cos(X(3))^2*cos(X(6))^3)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^3*sin(X(3))^2)/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^3*cos(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^3*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3))) + (m2*s2*sin(X(6))*(l2 - s2)*((v^2*cos(X(3))^2*cos(X(6))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (v^2*cos(X(6))^2*sin(X(3))^2*sin(X(6)))/(l1*(cos(X(3))*sin(X(6)) - cos(X(6))*sin(X(3)))^2) + (l2*v^2*cos(X(3))^2*cos(X(6))^2*sin(X(3)))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2 + (l2*v^2*cos(X(3))^2*sin(X(3))*sin(X(6))^2)/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))^2))/(l2*cos(X(3))*sin(X(6)) - l2*cos(X(6))*sin(X(3)))];
 
Q=Qg+Qf;
Cq=[1, 0,         s1*sin(X(3)),  0,  0,                    0;
    0, 1,        -s1*cos(X(3)),  0,  0,                    0;
    1, 0, -sin(X(3))*(l1 - s1), -1,  0,        -s2*sin(X(6));
    0, 1,  cos(X(3))*(l1 - s1),  0, -1,         s2*cos(X(6));
    0, 0,                    0,  1,  0, -sin(X(6))*(l2 - s2);
    0, 0,                    0,  0,  1,  cos(X(6))*(l2 - s2)];
A=[M Cq';
    Cq zeros(6,6)];

alph=50;
beta=50;
C=[X(1)-s1*cos(X(3));
    X(2)-s1*sin(X(3));
    X(1)+(l1-s1)*cos(X(3))-X(4)+s2*cos(X(6));
    X(2)+(l1-s1)*sin(X(3))-X(5)+s2*sin(X(6));
    X(4)+(l2-s2)*cos(X(6))-(v*t+cx);
    X(5)+(l2-s2)*sin(X(6))-cy];
C1=[X(7)+s1*X(9)*sin(X(3));
    X(8)-s1*X(9)*cos(X(3));
    X(7)-(l1-s1)*X(9)*sin(X(3))-X(10)-s2*X(12)*sin(X(6));
    X(8)+(l1-s1)*X(9)*cos(X(3))-X(11)-s2*X(12)*cos(X(6));
    X(10)-(l2-s2)*X(12)*sin(X(6))-v;
    X(11)+(l2-s2)*X(12)*cos(X(6))];
Gm=[s1*X(9)^2*cos(X(3));
    s1*X(9)^2*sin(X(3));
    -(l1-s1)*X(9)^2*cos(X(3))-s2*X(12)^2*cos(X(6));
    -(l1-s1)*X(9)^2*sin(X(3))-s2*X(12)^2*sin(X(6));
    -(l2-s2)*X(12)^2*cos(X(6));
    -(l2-s2)*X(12)^2*sin(X(6))]-2*alph*C1-(beta^2)*C;
RHS=[Q;
    Gm];
ACC=A\RHS;
dXdt=zeros(12,1);
dXdt(1:6)=X(7:12);
dXdt(7:12)=ACC(1:6);
    