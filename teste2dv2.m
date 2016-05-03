clear
close all
clc

load('paredes.mat');
IE = 1000;
JE = 1000;

olay = 302;

ic = 240+olay;
jc = 350+olay;

D = 0.03;   %Tamanho da c�lula

%Coordenadas dos receptores
Rx3x = 210+olay;
Rx3y = 300+olay;

Rx2x = 270+olay;
Rx2y = 300+olay;

c = 2.99792458e8;       %velocidade da luz
mi0 = ones(IE,JE)*4*pi*1e-7;          %permeabilidade do vacuo
eps0 = 1e-9/(36*pi)*epsilon;    %permissividade do vacuo
%mi0 = 1./(epsilon*c^2);

dt = 0.7*(D/(sqrt(2)*c));
nsteps = 500;
%t = dt*1:dt:nsteps*dt;

%fun��o geradora de pulso
Ap = 1;
t0 = 20*dt;
tau = 6*dt;
t = [1:nsteps]*dt;
g = -Ap*sqrt(2*exp(1)/tau^2)*(t-t0).*exp(-((t-t0)/tau).^2);

Ez = zeros(IE,JE);
Hx = zeros(IE,JE);
Hy = zeros(IE,JE);

figure(1)
plot(t,g)
title('Pulse in time')
xlabel('t')
ylabel('g(t)')

EzRx1 = zeros(1,nsteps);
EzRx2 = zeros(1,nsteps);
EzRx3 = zeros(1,nsteps);

tic
for k = 2:nsteps
    for j = 2:JE-1
        for i = 2:IE-1
            Ez(i,j) = Ez(i,j) + (dt/eps0(i,j))*(Hy(i,j)-Hy(i-1,j)-Hx(i,j)+Hx(i,j-1))/D;
        end
    end
    
    EzRx1(k) = Ez(ic,jc);
    EzRx2(k) = Ez(Rx2y,Rx2x);
    EzRx3(k) = Ez(Rx3y,Rx3x);
    
%     g = exp(-0.5*((t0-k)/tau)^2);
    Ez(ic,jc) = g(k);
    
    for j = 1:JE-1
        for i = 1:IE-1
            Hx(i,j) = Hx(i,j) + (dt/mi0(i,j))*(Ez(i,j)-Ez(i,j+1))/D;
            Hy(i,j) = Hy(i,j) + (dt/mi0(i,j))*(Ez(i+1,j)-Ez(i,j))/D;
        end
    end
    
    
    
    if k == 100
        salvoEz=Ez;
    end
    if k == 250
        salvoEz2=Ez;
    end
    if k == 350
        salvoEz3=Ez;
    end
    clc
    k
end
tempo = toc


figure(2)
subplot(2,2,1)
imagesc(abs(salvoEz(303:699,303:699)));
set(gca,'YDir','normal')

subplot(2,2,2)
imagesc(abs(salvoEz2(303:699,303:699)));
set(gca,'YDir','normal')

subplot(2,2,3)
imagesc(abs(salvoEz3(303:699,303:699)));
set(gca,'YDir','normal')

subplot(2,2,4)
imagesc(abs(Ez(303:699,303:699)))
set(gca,'YDir','normal')

figure(3)
subplot(3,1,1)
plot(t,EzRx1)
subplot(3,1,2)
plot(t,EzRx2)
subplot(3,1,3)
plot(t,EzRx3)

save('var2.mat', 'salvoEz', 'salvoEz2', 'salvoEz3', 'Ez', 'EzRx1', ...
    'EzRx2', 'EzRx3');


