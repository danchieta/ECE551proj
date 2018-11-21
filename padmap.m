function [newmap] = padmap(map, value, borderSize)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

newmap = ones(size(map)+borderSize*2)*value;
[h,w] = size(map);
newmap(borderSize+1: h+borderSize,borderSize+1: w+borderSize) = map;

end

