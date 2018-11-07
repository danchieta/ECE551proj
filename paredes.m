clear
close all
clc

im = imread('scene.png');

e1 = 5;     %permissividade das paredes externas
e2 = 4.2;   %permissividade das paredes internas
epsilon = ones(1000,1000); % inicia matriz com os valores de permissividade
%seguindo as orientaçõs do prof., faremos uma matriz maior do que o
%ambiente a ser estudado.

%O ambiente é quadrangular com aproximadamente 397 células de lado, então fiz
%uma matriz 1000 por 1000.

% Para converter as coordenada do mapa do prof. Muller para as do nosso,
% basta somar 302 aos valores de lá.

%PONHA SEUS LOOPS AQUI

for i=1:397
    for j=1:395
        if im(i,j)==1  %Parede externa
            epsilon(302+398-i,302+j) = e1;
        elseif im(i,j)==0  %parede interna
            epsilon(302+398-i,302+j) = e2;
        end
    end
end

figure
imagesc(1-epsilon); 
%imagesc(epsilon)
colormap(gray)
set(gca,'YDir','normal')

save('paredes.mat', 'epsilon');