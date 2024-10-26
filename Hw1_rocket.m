clc;clear;close all;

Gt = 0; #1 pentru traiectorie
Gd = 1; #1 pentru simulare
roVar = 1; #1 pentru densitatea aerului variabila cu altitudinea

g =  9.81;
D = 0.18;
m0 = 194 ; # masa rachetei
mc = 0.72 * m0; #masa combustibil
mr = m0 - mc;
mi = m0;
tao = 57;
a = 1; #daca se depaseste tao
q = mc / tao; #debit ardere
u = 3880; #viteza gazelor
rad = 0;
R = 8.3144 ;
E = 2.718;
Hn = 10.4;

v0 = 16;
alpha0 = 53;
eta = 1.81 * 1e-5;
b1 = 6.54 * eta * D;
c = 0.64;
ro0 = 1.22;
b2 = c * ro0 * D^2;
rad = 1;

t0 = 0;
vtao = 0;
# pentru aproximarea timpului max de zbor, presupunem masa rachetei este minimala
h0 = v0 * tao * cosd(alpha0) +((u * q * cosd(alpha0) / mr)-g) * tao^2/2;# h cand combustibilul se consuma
V = v0*cosd(alpha0) +(u * q * tao * cosd(alpha0)  / mr)- g * tao; #viteza cand combustibilul se termina
tz = 2 * V / g;
# dupa arderea combustibilului si deplasarea doar sub actiunea greutatii, pana la aceeasi inaltime de unde sa
#terminat combustibilul
t1 =(-V+ sqrt(V^2 + 2 * g *abs(h0))) /g;#restul de timp
tf = tao + tz + t1;
N =2000; #numarul de itervale pentru discretizarea timpului
t = linspace(t0,tf,N);
dt = t(2) - t(1);

vx = zeros(1,N);vy = vx;
x = zeros(1,N); y = x;
T = zeros(1,N); T(1) = 20;
P = zeros(1,N);P(1) = 101325;
vx(1) = v0 * sind(alpha0);
vy(1) = v0 * cosd(alpha0);

for i = 1:N-1
  if i * dt <= tao, mo = mr + mc - a * i *dt * q;
    else, a = 0;  vtao = vx(1)^2 + vy(i)^2;
  endif
  T(i+1) = T(i) - 6.5 * 1e-3 * y(i);
  P(i+1) = P(1) * E ^ (-1 * y(i) / Hn); # aproximarea presiunii in functie de altitudine: https://en.wikipedia.org/wiki/Density_of_air
  ro = 1e-3 * P(i) /(T(i) * R);
 b2 = c * ( (1 - roVar) * ro0 + roVar * ro) * D^2;

  rad = sqrt(vx(i)^2 + vy(i)^2);
  aux = 1 - dt * (b1 + b2 * rad - a * q -  a * q * u / rad)/m0;
  vx(i+1) = vx(i) * aux;
  vy(i+1) =  vy(i) * aux - g * dt;
  x(i+1) = x(i) + vx(i) *dt;
  y(i+1) = y(i) + vy(i) *dt;
    if y(i+1)<0, break; end
end

t = t(1:i) ; vx = vx(1:i); vy = vy(1:i); x = x(1:i); y = y(1:i);

tf = t(i);
b = x(i);
h = max(y);
tu = t (y==h);
tc = tf -tu;
mm = (mi + m0) / 2;
if tf < tao, Q = 1/2 *mm * (v0^2 - vx(i)^2 - vy(i)^2)  +  tf * q * (u^2);
  else, Q = 1/2 *mm * (v0^2 - vtao) + 1/2 * m0 * (vtao - vx(i)^2 - vy(i)^2)  + tao * q * (u^2);
endif
afis = ['Timpul de zbor: ', num2str(tf), ' s'];
disp(afis);
afis = ['Bataia proiectilului: ', num2str(b/1e3), ' km'];
disp(afis);
afis = ['Altitudinea max: ' , num2str(h/1e3), ' km'];
disp(afis);
afis = ['Timp urcare: ' , num2str(tu), ' s'];
disp(afis);
afis = ['Timp coborare: ' , num2str(tc), ' s'];
disp(afis);
afis = ['Energie disipata: ' , num2str(Q/1e6), ' MJ'];
disp(afis);

if Gt == 1 # forma traiectoriei
  figure(3);
  plot(x/1e3, y/1e3, 'k');
  xlabel('x(km)'); ylabel('y(km)');
  title('Curba balistica');
  axis equal; axis tight;
end

if Gd == 1 # dinamica in timp real a miscarii rachetei
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






