clear
close all
clc

table = [0 1 2;
    4.2 5 1];
epsilon = permittivitymap('scene.png',table);

% width of border arround the map that won't be shown but will exist for
% simulation purposes
padding = 302;

% map of relative permittivity
epsilon = padmap(epsilon,1,padding);

% dimensions of the simulated environment
env_size = size(epsilon);
[IE, JE] = size(epsilon);

% set here whether a video will be generated
video = 1;

% Coordinates of the source
source = sub2ind(size(epsilon),164+padding,295+padding);

D = 0.03;   % length of side of each cell

%Coordinates of receptors:
%Receptor 1
Rx3x = 210+padding;
Rx3y = 300+padding;

%Receptor 2
Rx2x = 270+padding;
Rx2y = 300+padding;

% Environment variables
c = 2.99792458e8;               % light speed
mu0 = 4*pi*1e-7;                % free space permeability
eps0 = 1/(c^2*mu0);
mir = ones(env_size)*4*pi*1e-7;    %permeabilidade do vacuo
epsr = 1e-9/(36*pi)*epsilon;    %permissividade do vacuo
%mi0 = 1./(epsilon*c^2);

dt = 0.7*(D/(sqrt(2)*c));
nsteps = 800;
%t = dt*1:dt:nsteps*dt;

%fun��o geradora de pulso
Ap = 1;                 %Amplitude m�xima
t0 = 20*dt;             %Tempo onde o pulso est� localizado
tau = 6*dt;             %Largura do pulso
t = [1:nsteps]*dt;      %vetor de tempo
g = -Ap*sqrt(2*exp(1)/tau^2)*(t-t0).*exp(-((t-t0)/tau).^2); %fun��o

%iniciando vetores antes do loop
Ez = zeros(env_size);
Hx = zeros(env_size);
Hy = zeros(env_size);

figure(1)
%Plota a fun��o geradora de pulso
plot(t/1e-9,g)
title('Pulse in time')
xlabel('t (ns)')
ylabel('g(t)')

if video
    %Objeto para exportar v�deo com os quadros
    v = VideoWriter('FDTD.avi');
    v.Quality = 100;
    v.FrameRate = 30;
    open(v);
end

%iniciando vetores antes do loop
EzRx1 = zeros(1,nsteps);
EzRx2 = zeros(1,nsteps);
EzRx3 = zeros(1,nsteps);


fcounter = 1;

sez = zeros(nsteps,2);
tic
% loop through time
for k = 2:nsteps
    clc
    k
    sez(k,:) = size(Ez);
    Ez(2:end,2:end) = Ez(2:end,2:end) ...
        + (dt./epsr(2:end,2:end))...
        .*(Hy(2:end,2:end)-Hy(1:end-1,2:end)...
        -Hx(2:end,2:end)+Hx(2:end,1:end-1))/D;
    
    Ez(source) = g(k);   % Excita��o na fonte do pulso
    
    EzRx1(k) = Ez(source);       % Campo detectado na fonte
    EzRx2(k) = Ez(Rx2y,Rx2x);   % Campo detectado em Rx2
    EzRx3(k) = Ez(Rx3y,Rx3x);   % Campo detectado em Rx3
    
    Hx(:, 1:end-1) = Hx(:, 1:end-1) - (dt./mir(:, 1:end-1))...
        .*(Ez(:, 2:end) - Ez(:, 1:end-1))/D;
    Hy(1:end-1, :) = Hy(1:end-1, :) + (dt./mir(1:end-1, :))...
        .*(Ez(2:end, :) - Ez(1:end-1, :))/D;
    
   
    if (k==2 || ~(rem(k,5)))&& video
        % capturing video frames
        figure(5)
        imagesc(abs(Ez(padding+1:end-padding,padding+1:end-padding)))
        axis('equal')
        title(['time: ' num2str(k*dt/1e-9) 'ns'])
        
        set(gcf,'Position', [100 100 720 480])

        M(fcounter) = getframe(gcf);
        fcounter  = fcounter +1;
        %close(gcf)
        
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
end
tempo = toc


writeVideo(v,M);
close(v); % Grava o v�deo  no arquivo

figure(2)
subplot(2,2,1)
imagesc(abs(salvoEz(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Passo 100; ' num2str(100*dt/1e-9) 'ns'])
axis('equal')

subplot(2,2,2)
imagesc(abs(salvoEz2(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Passo 250; ' num2str(250*dt/1e-9) 'ns'])
axis('equal')

subplot(2,2,3)
imagesc(abs(salvoEz3(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Passo 350; ' num2str(350*dt/1e-9) 'ns'])
axis('equal')

subplot(2,2,4)
imagesc(abs(Ez(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Passo 500; ' num2str(500*dt/1e-9) 'ns'])
axis('equal')

figure(3)
subplot(3,1,1)
plot(t/1e-9,EzRx1)
title('Intensidade de campo el�trico nos receptores')
xlabel('t (ns) - receptor 1')
ylabel('Ez')
subplot(3,1,2)
plot(t,EzRx2)
xlabel('t (ns) - receptor 2')
ylabel('Ez')
subplot(3,1,3)
plot(t,EzRx3)
xlabel('t (ns) - receptor 3')
ylabel('Ez')

figure(4)
imagesc(abs(salvoEz3(padding+1:end-padding,padding+1:end-padding)));
colorbar('EastOutside')
title(['Passo 350; ' num2str(350*dt/1e-9) 'ns'])
axis('equal')


save('var2.mat', 'salvoEz', 'salvoEz2', 'salvoEz3', 'Ez', 'EzRx1', ...
    'EzRx2', 'EzRx3');


