clc;clear;close all;

Gv = 1;
Gc = 1;
Gt = 1;
Gd = 1;
if Gd == 1, Gv = 0, Gc = 0, Gt = 0; end

g =  9.80665;
ro = 7850;
r = 0.13;
m =4/3 * pi * r^3 * ro;
G = m * g;

v0 = 1100;
alpha0 = 43;
eta = 1.81 * 1e-5;
b1 = 6 * pi * eta * r;
c = 0.469;
ro0 = 1.22;
b2 = c * 4 * pi * r^2 *ro0 /2;

t0 = 0;
tf = 2 * v0 /g*sind(alpha0);
N = 1500;
t = linspace(t0,tf,N);
dt = t(2) - t(1);

vx = zeros(1,N);vy = vx;
x = zeros(1,N); y = x;
vx(1) = v0 * sind(alpha0);
vy(1) = v0 * cosd(alpha0);

for i = 1:N-1
  aux = 1 - dt * (b1 + b2 * sqrt(vx(i)^2 + vy(i)^2))/m;
  vx(i+1) = vx(i) * aux;
  vy(i+1) =  vy(i) * aux - g * dt;
  x(i+1) = x(i) + vx(i) *dt;
  y(i+1) = y(i) + vy(i) *dt;
    if y(i+1)<0, break; end
end

t = t(1:i) ; vx = vx(1:i); vy = vy(1:i); x = x(1:i); y = y(1:i);

if Gv == 1
  figure(1);
  plot(t, vx, '-r', t, vy, '-b');
  xlabel('t(s)') ;ylable = ('v(m/s)'); grid;
  title('Vitezele ca functie de timp');
  legend('vx', 'vy');
end

if Gc == 1
   figure(2);
   plot(t, x/ 1e3, '-r', t, y/1e3, '-b');
   xlabel('(t(s)'); ylabel('coord(km)');grid;
   title('Coordonatele ca functie de timp');
   legend('x','y');
end

if Gt == 1
  figure(3);
  plot(x/1e3, y/1e3, '-k', 'LineWith', 2);
  xlabel('x(km)'); ylabel('y(km)');
  title('Curba balistica');
  axis equal; axis tight;
end

tf = t(i);
b = x(i);
h = max(y);
tu = t (y==h);
tc = tf -tu;
Q = 1/2 * m * (v0^2 - vx(i)^2 - vy(i)^2);
afis = ['Timpul de zbor :', num2str(tf), ' s'];
disp(afis);
afis = ['Bataia proiectilului :', num2str(b/1e3), ' km'];
disp(afis);
afis = ['Altitudinea max:' , num2str(h/1e3), ' km'];
disp(afis);



if Gd == 1
  figure(4);
  set(4, 'Position', [50 50 850 600]);
  tic; simt=0;
 while simt < tf
   plot(x/1e3,y/1e3, '-c');hold on;
   xlabel('x(km)'); ylabel('y(km)');
   grid;
   title('Simularea miscarii');
   axis equal; axis tight;
   index = abs(t - simt) == min(abs(t-simt));
   plot(x(index)/1e3, y(index)/1e3, '.b', 'MarkerSize', 10);
   hold off

   text(b/2/1e3, h/3/1e3, ['vx=', num2str(round(vx(index))), ' m/s']);
   text((b/2 - b/5)/1e3 , h/3/1e3,['t=', num2str(round(t(index))), ' s']);
   text((b/2+b/5)/1e3, h/3/1e3, ['vy=', num2str(round(vy(index))), ' m/s']);
   pause(1e-3);
   simt = toc;
 end
end






