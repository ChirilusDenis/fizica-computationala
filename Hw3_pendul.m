clc;clear; close all;

L = 1;  % lungime pendul
Lm = L * 1.1;  % coordonate maxime pentru axele simularii
g = 9.80665;
m = 0.9;  % masa obiectului
omega = sqrt(g/L);  % pulsatia
T = 2 * pi / omega;  % calcularea unei perioade
theta0 = 45 * pi / 180;  % unghiul initial in radiani
vi = 120 * pi / 180;  % viteza initilala in radiani pe secunda

ti=0; tf=5*T;  % afisarea a 5 perioade
N=100000;  % numar diviziuni
t=linspace(ti,tf,N);  % spatiul timpului
dt=t(2)-t(1);  % diferenta intre 2 momente de timp

theta = zeros(1, N);  % vectorul unghiurilor
vu = theta;  % vectorul vitezelor initiale
theta(1) = theta0;
theta(2) = theta0 + vi * dt;
vu(1) = vi;

for i = 2:N-1
  % formula finala rezultata din lagrangian
  % theta'' + g/L * sin(theta) = 0;
  vu(i) = (theta(i) - theta(i-1)) / dt;
  e = - g / L * sin(theta(i));
  theta(i + 1) = 2 * theta(i) - theta(i - 1)+dt^2 * e;
endfor

x = L * sin(theta);  % coordonata x
y =  - L * cos(theta);  % coordonata y
T = m * vu .^2 * L^2  / 2;  % energia cinetica
U = L * m * g * (1 - cos(theta));  % energia potentiala cu pendul legat la h = L de axa x = 0
H = T + U;  % hamiltonian, energia totala

tic; simt = 0;
while simt < tf
   j=abs(t-simt)==min(abs(t-simt));  % moment timp cel mai apropiat de timpul de simulare
   plot([0 x(j)],[0 y(j)],'-k','LineWidth',1); hold on;  % tija pendulului
   axis square;
   plot(x(j),y(j),'.b','MarkerSize',30);  % obiectul pendulului
   % text(3/5 * Lm, 3/5 * Lm,['E = ',num2str(round(H(j))),' J']);
   text(3/5 * Lm, 5/6 * Lm,['E = ',num2str(H(j)),' J']);  % afisarea energiei
   axis([-Lm , Lm, -Lm, Lm]); hold off;
   simt = toc;
   pause(1e-6);
endwhile
