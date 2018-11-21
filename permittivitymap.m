function [map] = permittivitymap(img_name, epsilon_table)
%permittivitymap This function creates the map of electrict permitivity of
%the environment to be simulated.
%   img_name is the name of the image on which the map will be based.
%   epsilon_table 2-row matrix relating the values of the image on first
%   row with dieletric constant in the second row.

img = imread(img_name); %read the image
map = zeros(size(img)); %create the map

[~,N] = size(epsilon_table);

% loop through each column in table
for i=1:N 
    % relate the value image in with electric permitivity.
    epsilon = (img==epsilon_table(1,i))*epsilon_table(2,i);
    map = map + epsilon;
    
end

end

