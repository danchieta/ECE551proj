IE = int16(300);
JE = int16(300);

ic = IE/2;
jc = JE/2;

D = 0.015;

c = 2.99792458e8;       %velocidade da luz
mi0 = pi*4e-7;          %condutividade do vacuo
eps0 = 1e-9/(36*pi);    %permissividade do vacuo

dt = 0.99*(D/(sqrt(2)*c));
nsteps = 500;
%t = dt*1:dt:nsteps*dt;

%fun��o geradora de pulso
A0 = 5;
t0 = 20;
tau = 6;
%g = A0*exp(-(t-t0).^2/tau^2);
%g = zeros(length(t),1);

Ez = zeros(IE,JE);
Hx = zeros(IE,JE);
Hy = zeros(IE,JE);
% figure(1)
% plot(k,g)
% title('Pulse in time')
% xlabel('t')
% ylabel('g(t)')

tic
for k = 2:nsteps
    for j = 2:JE-1
        for i = 2:IE-1
            Ez(i,j) = Ez(i,j) + (dt/eps0)*(Hy(i,j)-Hy(i-1,j)-Hx(i,j)+Hx(i,j-1))/D;
        end
    end
    
    g = exp(-0.5*((t0-k)/tau)^2);
    Ez(ic,jc) = g;
    
    for j = 1:JE-1
        for i = 1:IE-1
            Hx(i,j) = Hx(i,j) + (dt/mi0)*(Ez(i,j)-Ez(i,j+1))/D;
            Hy(i,j) = Hy(i,j) + (dt/mi0)*(Ez(i+1,j)-Ez(i,j))/D;
        end
    end
    
    if k == 125
        salvoEz=Ez;
    end
    if k == 250
        salvoEz2=Ez;
    end
    if k == 375
        salvoEz3=Ez;
    end
    clc
    progress = 100*k/nsteps
end
tempo = toc


figure(2)
subplot(2,2,1)
imagesc(abs(salvoEz));
axis('equal')

subplot(2,2,2)
imagesc(abs(salvoEz2));

subplot(2,2,3)
imagesc(abs(salvoEz3));

subplot(2,2,4)
imagesc(abs(Ez))




