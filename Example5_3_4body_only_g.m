function Example5_3
clear all
global m1 m2 m3 m4 I1 I2 I3 I4 g l1 l2 l3 l4 s1 s2 s3 s4 v  cx cy u Dx Dy
m1=1;
m2=1;
m3=1;
m4=1;
l1=1;
l2=1;
l3=1;
l4=1;
s1=0.5;
s2=0.5;
s3=0.5;
s4=0.5;
g=9.81;
I1=m1*l1^2/12;
I2=m2*l2^2/12;
I3=m3*l3^2/12;
I4=m4*l4^2/12;
v=0.5;
cx=l1*cos(pi/4)+l2*cos(pi*2/3);
cy=l1*sin(pi/4)+l2*sin(pi*2/3);
u=0.5;
Dx=0;
Dy=0;

TSPAN=[0:1E-2:2];
X0=[s1*cos(pi/4) s1*sin(pi/4) pi/4 l1*cos(pi/4)+s2*cos(pi*2/3) l1*sin(pi/4)+s2*sin(pi*2/3) pi*2/3 l4*cos(pi/4)+s3*cos(pi*2/3) l4*sin(pi/4)+s3*sin(pi*2/3) pi*2/3 s4*cos(pi/4) s4*sin(pi/4) pi/4 0 0 0 0 0 0 0 0 0 0 0 0];
[t,X]=ode45(@func,TSPAN,X0);

for i=1:length(X)
    figure(1);
    plot([0,0],[0,0]);
    grid;
    axis square;
    xmin=-3;xmax=3;ymin=-3;ymax=3;
    axis([xmin,xmax,ymin,ymax]);
    hold on;
    xQp=X(i,1)-s1*cos(X(i,3));
    yQp=X(i,2)-s1*sin(X(i,3));
    xPp1=X(i,1)+(l1-s1)*cos(X(i,3));
    yPp1=X(i,2)+(l1-s1)*sin(X(i,3));
    xPp2=X(i,4)-s2*cos(X(i,6));
    yPp2=X(i,5)-s2*sin(X(i,6));
    xRp=X(i,4)+(l2-s2)*cos(X(i,6));
    yRp=X(i,5)+(l2-s2)*sin(X(i,6));
    xRs=X(i,7)+(l3-s3)*cos(X(i,9));
    yRs=X(i,8)+(l3-s3)*sin(X(i,9));
    xPs2=X(i,7)-s3*cos(X(i,9));
    yPs2=X(i,8)-s3*sin(X(i,9));
    xPs1=X(i,10)+(l4-s4)*cos(X(i,12));
    yPs1=X(i,11)+(l4-s4)*sin(X(i,12));
    xQs=X(i,10)-s4*cos(X(i,12));
    yQs=X(i,11)-s4*sin(X(i,12));
    plot([xQp,xPp1],[yQp,yPp1],'linewidth',4);
    plot([xRp,xPp2],[yRp,yPp2],'linewidth',4);
    plot([xQs,xPs1],[yQs,yPs1],'linewidth',4);
    plot([xRs,xPs2],[yRs,yPs2],'linewidth',4);
    
    drawnow
    hold off
    set(1,'doublebuffer','on')
end

figure(2);
plot(t,X(:,5)+(l2-s2)*sin(X(:,6)))
xlabel('time','FontSize',18);
ylabel('Ry','FontSize',18);

function dXdt=func(t,X)
global m1 m2 m3 m4 I1 I2 I3 I4 g l1 l2 l3 l4 s1 s2 s3 s4 v  cx cy u Dx Dy
M1=[m1 0 0;
    0 m1 0;
    0 0 I1];
M2=[m2 0 0;
    0 m2 0;
    0 0 I2];
M3=[m3 0 0;
    0 m3 0;
    0 0 I3];
M4=[m4 0 0;
    0 m4 0;
    0 0 I4];
Q=[0;
    -m1*g;
    0;
    0;
    -m2*g;
    0;
    0;
    -m3*g;
    0;
    0;
    -m4*g;
    0];

Cq1=[1 0 s1*sin(X(3));
    0 1 -s1*cos(X(3));
    1 0 -(l1-s1)*sin(X(3));
    0 1 (l1-s1)*cos(X(3));
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0];
Cq2=[0 0 0;
    0 0 0;
    -1 0 -s2*sin(X(6));
    0 -1 s2*cos(X(6));
    1 0 -(l2-s2)*sin(X(6));
    0 1 (l2-s2)*cos(X(6));
    0 0 0;
    0 0 0];
Cq3=[0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
    -1 0 (l3-s3)*sin(X(9));
    0 -1 -(l3-s3)*cos(X(9));
    1 0 s3*sin(X(9));
    0 1 -s3*cos(X(9))];
Cq4=[0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    -1 0 (l4-s4)*sin(X(12));
    0 -1 -(l4-s4)*cos(X(12))]; 
A=[M1 zeros(3,3) zeros(3,3) zeros(3,3)   Cq1';
    zeros(3,3) M2 zeros(3,3) zeros(3,3)  Cq2';
    zeros(3,3) zeros(3,3) M3 zeros(3,3)  Cq3';
    zeros(3,3) zeros(3,3) zeros(3,3) M4  Cq4';
    Cq1        Cq2        Cq3       Cq4  zeros(8,8)];

alph=50;
beta=50;
C=[X(1)-s1*cos(X(3));
    X(2)-s1*sin(X(3));
    X(1)+(l1-s1)*cos(X(3))-X(4)+s2*cos(X(6));
    X(2)+(l1-s1)*sin(X(3))-X(5)+s2*sin(X(6));
    X(4)+(l2-s2)*cos(X(6))-X(7)-(l3-s3)*cos(X(9));
    X(5)+(l2-s2)*sin(X(6))-X(8)-(l3-s3)*sin(X(9));
    X(7)-s3*cos(X(9))-X(10)-(l4-s4)*cos(X(12));
    X(8)-s3*sin(X(9))-X(11)-(l4-s4)*sin(X(12))];
    
C1=[X(13)+s1*sin(X(3))*X(15);
    X(14)-s1*cos(X(3))*X(15);
    X(13)-(l1-s1)*sin(X(3))*X(15)-X(16)-s2*sin(X(6))*X(18);
    X(14)+(l1-s1)*cos(X(3))*X(15)-X(17)+s2*cos(X(6))*X(18);
    X(16)-(l2-s2)*sin(X(6))*X(18)-X(19)+(l3-s3)*sin(X(9))*X(21);
    X(17)+(l2-s2)*cos(X(6))*X(18)-X(20)-(l3-s3)*cos(X(9))*X(21);
    X(19)+s3*sin(X(9))*X(21)-X(22)+(l4-s4)*sin(X(12))*X(24);
    X(20)-s3*cos(X(9))*X(21)-X(23)-(l4-s4)*cos(X(12))*X(24)];    
    

Gm=[s1*sin(X(3))*X(15)^2;
    s1*cos(X(3))*X(15)^2;
    -(l1-s1)*cos(X(3))*X(15)^2-s2*cos(X(6))*X(18)^2;
    -(l1-s1)*sin(X(3))*X(15)^2-s2*sin(X(6))*X(18)^2;
    -(l2-s2)*cos(X(6))*X(18)^2+(l3-s3)*cos(X(9))*X(21)^2;
    -(l2-s2)*sin(X(6))*X(18)^2+(l3-s3)*sin(X(9))*X(21)^2;
    s3*cos(X(9))*X(21)^2+(l4-s4)*cos(X(12))*X(24)^2;
    s3*sin(X(9))*X(21)^2+(l4-s4)*sin(X(12))*X(24)^2]-2*alph*C1-(beta^2)*C;
RHS=[Q;
    Gm];
ACC=A\RHS;
dXdt=zeros(24,1);
dXdt(1:12)=X(13:24);
dXdt(13:24)=ACC(1:12);

    