clc; clear all; close all;

xmax = 1;
ymax = 1;
vmin = 0;
vmax = 2000;

kB=1.38e-23;
u=1.66e-27;
m0 = 32 * u;  % masa unei particule
T = 300;  % temperatura
nrp =100;  % numar puncte/ particule
count = 1;
nr = 0;
v1 = 0;
v2 = 100;
speeds = zeros(1, 1000);

while v2 <= vmax
nr = round(integral(@(v) (4*pi*v.^2*(m0/(2*pi*kB*T))^(3/2).*exp(-m0*v.^2/(2*kB*T))), v1, v2) * 1000);
% nr este frecventa relativa a vitezelor intre v1 si v2 in distributia Maxwell
if nr != 0
speeds(count : (count + nr)) = (v1 + v2) / 2;
endif
count = count + nr;  % de unde trebuie continuarea de completare a vectorului de speeds
v1 = v2;
v2 = v2 +100;  % urmatorul interval de viteze
endwhile

x = rand(1, nrp);
y = rand(1, nrp);
x = x * xmax;
y = y * ymax;
% ^ punere de puncte random in patratul [0, xmax]x[0, ymax]
alpha = rand(1, nrp);
alpha = alpha .*2 .*pi;
% ^ unghiuri random ale vitezelor
v = zeros(1, nrp);
v = speeds(round(rand(1, 100)*1000));
% ^ viteze random din distibutie
vx = v.*cos(alpha);
vy = v.*sin(alpha);

figure

for i = 1:100
  x =  x + vx / 5e4;
  y =  y + vy / 5e4;
% ^ urmatoare pozitie pentru particule cu viteze diminuate pentru vizualizare

  p = find(x >= 1);
  vx(p) = vx(p) .* (-1);
  p = find(x <= 0);
  vx(p) = vx(p) .* (-1);
  p = find(y >= 1);
  vy(p) = vy(p) .* (-1);
  p = find(y <= 0);
  vy(p) = vy(p) .* (-1);
% ^ particulele care se lovesc de pereti, ricoseaza inapoi

  plot(x, y, '.', "Markersize", 10);
  axis([0, xmax, 0, ymax], "manual");
  pause(0.1);
endfor

