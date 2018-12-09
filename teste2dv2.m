clear
% close all
clc

% filename of the scene image file
% use a grayscale png file preferably
scene = 'scene2.png';
table = [0 255;
     5 1];
epsilon = permittivitymap(scene,table);


% width of border arround the map that won't be shown but will exist for
% simulation purposes
padding = 302;

% map of relative permittivity
epsilon = padmap(epsilon,1,padding);

% dimensions of the simulated environment
env_size = size(epsilon);
[IE, JE] = size(epsilon);

% set here whether a video will be generated or not
video = 1;

% Coordinates of the source
Tx = 35;
Ty = 166;
source = sub2ind(size(epsilon),Tx+padding,Ty+padding);

D = 0.03;   % length of side of each cell

%Coordinates of receptors:
%Receptor 1
Rx3x = 164+padding;
Rx3y = 161+padding;

%Receptor 2
Rx2x = 163+padding;
Rx2y = 32+padding;

% Environment variables
c = 2.99792458e8; % light speed
mu0 = 4*pi*1e-7; % free space permeability
eps0 = 1/(c^2*mu0); % free space permittivity            

% calculating actual permittivity and permeability from relative values
mir = ones(env_size)*4*pi*1e-7; 
epsr = eps0*epsilon;

dt = 0.7*(D/(sqrt(2)*c));
nsteps = 800;
%t = dt*1:dt:nsteps*dt;

% pulse function
t = [0:nsteps-1]*dt;
% wavelength is 10 times size of a cell
lambda = 10*D;
f = c/lambda;
% duration of pulse
cycles = 5;
duration = ceil((cycles/f)/dt);
% pulse itself
g = sin(2*pi*f*[0:duration-1]*dt).*hamming(duration)';
g = [g zeros(1,nsteps - length(g))];

% initializin arrays before loop
Ez = zeros(env_size);
Hx = zeros(env_size);
Hy = zeros(env_size);
EzRx1 = zeros(1,nsteps);
EzRx2 = zeros(1,nsteps);
EzRx3 = zeros(1,nsteps);


figure(1)
imshow(imread(scene))
hold on
plot(Ty, Tx, 'bx')
text(Ty-39, Tx, 'Tx1/Rx1')
plot(Rx3y - padding, Rx3x - padding, 'bo')
text(Rx3y - padding -5, Rx3x - padding - 10, 'Rx3')
plot(Rx2y - padding, Rx2x - padding, 'bo')
text(Rx2y - padding + 5, Rx2x - padding, 'Rx2')

% plot pulse generating function
figure(2)
plot(t/1e-9,g)
title('Pulse in time')
xlabel('t (ns)')
ylabel('g(t)')

if video
    % declare video object
    v = VideoWriter('FDTD.avi');
    v.Quality = 100;
    v.FrameRate = 30;
    open(v);
end


snapshots = round([0.25:.25:1]*nsteps);
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
    
    Ez(source) = g(k);   % excitation at source
    
    % capturing detected field at transmitter and receptors locations
    EzRx1(k) = Ez(source+1);   
    EzRx2(k) = Ez(Rx2x,Rx2y);
    EzRx3(k) = Ez(Rx3x,Rx3y);
    
    Hx(:, 1:end-1) = Hx(:, 1:end-1) - (dt./mir(:, 1:end-1))...
        .*(Ez(:, 2:end) - Ez(:, 1:end-1))/D;
    Hy(1:end-1, :) = Hy(1:end-1, :) + (dt./mir(1:end-1, :))...
        .*(Ez(2:end, :) - Ez(1:end-1, :))/D;
    
   
    if (k==2 || ~(rem(k,5)))&& video
        % capturing video frames
        figure(3)
        imagesc(abs(Ez(padding+1:end-padding,padding+1:end-padding)))
        axis('equal')
        title(['time: ' num2str(k*dt/1e-9) 'ns'])
        xlabel('y Cells')
        ylabel('x Cells')
        
        set(gcf,'Position', [100 100 720 480])

        M(fcounter) = getframe(gcf);
        fcounter  = fcounter +1;
        
    end
    
    % Save some time samples of Ez
    if k == snapshots(1)
        salvoEz=Ez;
    end
    if k == snapshots(2)
        salvoEz2=Ez;
    end
    if k == snapshots(3)
        salvoEz3=Ez;
    end
end
tempo = toc


writeVideo(v,M);
close(v); % writes video file

figure(4)
subplot(2,2,1)
imagesc(abs(salvoEz(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Step ' num2str(snapshots(1)) '; ' num2str(snapshots(1)*dt/1e-9) 'ns'])
axis('equal')
xlabel('y Cells')
ylabel('x Cells')

subplot(2,2,2)
imagesc(abs(salvoEz2(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Step ' num2str(snapshots(2)) '; ' num2str(snapshots(2)*dt/1e-9) 'ns'])
axis('equal')
xlabel('y Cells')
ylabel('x Cells')

subplot(2,2,3)
imagesc(abs(salvoEz3(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Step ' num2str(snapshots(3)) '; ' num2str(snapshots(3)*dt/1e-9) 'ns'])
axis('equal')
xlabel('y Cells')
ylabel('x Cells')

subplot(2,2,4)
imagesc(abs(Ez(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Step ' num2str(snapshots(4)) '; ' num2str(snapshots(4)*dt/1e-9) 'ns'])
axis('equal')
xlabel('y Cells')
ylabel('x Cells')

figure(5)
subplot(3,1,1)
plot(t/1e-9,EzRx1)
title('Receptor 1')
xlabel('t (ns)')
ylabel('Ez (V/m)')

subplot(3,1,2)
plot(t,EzRx2)
title('Receptor 2')
xlabel('t (ns)')
ylabel('Ez (V/m)')

subplot(3,1,3)
plot(t,EzRx3)
title('Receptor 3')
xlabel('t (ns)')
ylabel('Ez (V/m)')

figure(6)
imagesc(abs(salvoEz2(padding+1:end-padding,padding+1:end-padding)))
colorbar('EastOutside')
title(['Step ' num2str(snapshots(2)) '; ' num2str(snapshots(2)*dt/1e-9) 'ns'])
axis('equal')
xlabel('y Cells')
ylabel('x Cells')


save('var2.mat', 'salvoEz', 'salvoEz2', 'salvoEz3', 'Ez', 'EzRx1', ...
    'EzRx2', 'EzRx3');


