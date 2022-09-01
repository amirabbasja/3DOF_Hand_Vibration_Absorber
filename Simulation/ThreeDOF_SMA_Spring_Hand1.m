clc;
clear;

syms th1(t) th2(t) xd(t) T1 T2 F_SMA
m1=1; m2 =1;md =1;g =1;c1 =1;c2 =1;cd =1;kd =1;l1 =1;l2 =1;l3 =1;lc1=1; lc2 =1;lc3=1; ld = 1;Kt1 = 1; Kt2 = 1;

T1 = sin(t);
T2 = cos(t);
F_SMA = 0;

eq1 = g* cos(th1)*lc1*m1+ld*md*cos(th1)-sin(th1)*md*xd+ Kt1*th1+ld*md*diff(xd,2)+lc1^2*m1*diff(th1,2)+ ...
    l1^2*m2*diff(th1,2)+ld^2*md*diff(th1,2)+m2*l1*(g*cos(th1)+lc2*(sin(th1-th2)*diff(th2,1)^2+cos(th1-th2)*diff(th2,2)))...
     == T1 - c1 * diff(th1,1);

eq2 = Kt2*th2+ lc2*m2*(g*cos(th2)+l1*(-sin(th1-th2)*diff(th1,1)^2+cos(th1-th2)*diff(th1,2))) + lc1^2*m2*diff(th2,2)...
    == T2 - c2*diff(th2,1);

eq3 = md*(diff(xd,2)+ld*diff(th1,2)) + md*cos(th1) == F_SMA - cd*diff(xd,1);

[V,S] = odeToVectorField(eq1, eq2, eq3)

M = matlabFunction(V,'vars', {'t','Y'});

  xd_0 = 0;
 Dxd_0 = 0;
 th1_0 = 0;
Dth1_0 = 0;
 th2_0 = 0;
Dth2_0 = 0;

sol = ode45(M,[0 80],[xd_0 Dxd_0 th1_0 Dth1_0 th2_0 Dth2_0]);

x = linspace(0,80,100);
y = deval(sol,x,1);
plot(x,y)

