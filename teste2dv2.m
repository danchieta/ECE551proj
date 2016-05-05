clear
close all
clc

load('paredes.mat');

IE = 1000;
JE = 1000;

video = 1;

olay = 302;

ic = 240+olay;
jc = 350+olay;

D = 0.03;   %Tamanho da célula

%Coordenadas dos receptores:
%Receptor 1
Rx3x = 210+olay;
Rx3y = 300+olay;

%Receptor 2
Rx2x = 270+olay;
Rx2y = 300+olay;

%Variáveis do ambiente
c = 2.99792458e8;               %velocidade da luz
mi0 = ones(IE,JE)*4*pi*1e-7;    %permeabilidade do vacuo
eps0 = 1e-9/(36*pi)*epsilon;    %permissividade do vacuo
%mi0 = 1./(epsilon*c^2);

dt = 0.7*(D/(sqrt(2)*c));
nsteps = 500;
%t = dt*1:dt:nsteps*dt;

%função geradora de pulso
Ap = 1;                 %Amplitude máxima
t0 = 20*dt;             %Tempo onde o pulso está localizado
tau = 6*dt;             %Largura do pulso
t = [1:nsteps]*dt;      %vetor de tempo
g = -Ap*sqrt(2*exp(1)/tau^2)*(t-t0).*exp(-((t-t0)/tau).^2); %função

%iniciando vetores antes do loop
Ez = zeros(IE,JE);
Hx = zeros(IE,JE);
Hy = zeros(IE,JE);

figure(1)
%Plota a função geradora de pulso
plot(t,g)
title('Pulse in time')
xlabel('t')
ylabel('g(t)')

if video
    %Objeto para exportar vídeo com os quadros
    vidObj = VideoWriter('FDTD.avi');
    vidObj.Quality = 100;
    vidObj.FrameRate = 30;
    open(vidObj);
end

%iniciando vetores antes do loop
EzRx1 = zeros(1,nsteps);
EzRx2 = zeros(1,nsteps);
EzRx3 = zeros(1,nsteps);

tic
%loop no tempo
for k = 2:nsteps
    for j = 2:JE-1
        for i = 2:IE-1
            %Cálculo do campo E para todos os pontos no espaço
            Ez(i,j) = Ez(i,j) + (dt/eps0(i,j))*(Hy(i,j)-Hy(i-1,j)-Hx(i,j)+Hx(i,j-1))/D;
        end
    end
    
    Ez(ic,jc) = g(k);   % Excitação na fonte do pulso
    
    EzRx1(k) = Ez(ic,jc);       % Campo detectado na fonte
    EzRx2(k) = Ez(Rx2y,Rx2x);   % Campo detectado em Rx2
    EzRx3(k) = Ez(Rx3y,Rx3x);   % Campo detectado em Rx3
    
%     g = exp(-0.5*((t0-k)/tau)^2);

    for j = 1:JE-1
        for i = 1:IE-1
            % Cálculo do campo H para todos os pontos no espaço
            Hx(i,j) = Hx(i,j) + (dt/mi0(i,j))*(Ez(i,j)-Ez(i,j+1))/D;
            Hy(i,j) = Hy(i,j) + (dt/mi0(i,j))*(Ez(i+1,j)-Ez(i,j))/D;
        end
    end
    
    if (k==2 || ~(rem(k,5)))&& video
        % Aqui são feitos os quadros pro vídeo
        h= figure;
        imagesc(abs(Ez(303:699,303:699)))
        set(gca,'YDir','normal')
        title(['time: ' num2str(k*dt/1e-9) 'ns'])
        writeVideo(vidObj, getframe(h));
        close 
        
    end
    
    % Salva algumas amostras de Ez no tempo
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

close(vidObj); % Grava o vídeo  no arquivo

figure(2)
subplot(2,2,1)
imagesc(abs(salvoEz(303:699,303:699)));
colorbar('EastOutside')
title(['Passo 100; ' num2str(100*dt/1e-9) 'ns'])
set(gca,'YDir','normal')

subplot(2,2,2)
imagesc(abs(salvoEz2(303:699,303:699)));
colorbar('EastOutside')
title(['Passo 250; ' num2str(250*dt/1e-9) 'ns'])
set(gca,'YDir','normal')

subplot(2,2,3)
imagesc(abs(salvoEz3(303:699,303:699)));
colorbar('EastOutside')
title(['Passo 350; ' num2str(350*dt/1e-9) 'ns'])
set(gca,'YDir','normal')

subplot(2,2,4)
imagesc(abs(Ez(303:699,303:699)))
colorbar('EastOutside')
title(['Passo 500; ' num2str(500*dt/1e-9) 'ns'])
set(gca,'YDir','normal')

figure(3)
subplot(3,1,1)
plot(t/1e-9,EzRx1)
title('Intensidade dos campos elétricos nos receptores')
xlabel('t (ns) - receptor 1')
subplot(3,1,2)
plot(t,EzRx2)
xlabel('t (ns) - receptor 2')
subplot(3,1,3)
plot(t,EzRx3)
xlabel('t (ns) - receptor 3')

save('var2.mat', 'salvoEz', 'salvoEz2', 'salvoEz3', 'Ez', 'EzRx1', ...
    'EzRx2', 'EzRx3');


