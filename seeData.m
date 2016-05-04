close all
load('var1.mat')
im = imread('scene.png')

figure(2)
subplot(2,2,1)
imshow(im)
%hold on
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
plot(EzRx1)
subplot(3,1,2)
plot(EzRx2)
subplot(3,1,3)
plot(EzRx3)