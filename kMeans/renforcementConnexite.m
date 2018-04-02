function [ labels ] = renforcementConnexite(labels,Centers)
%CONNEXITE2 Summary of this function goes here
%   Detailed explanation goes here

for k = 1:size(Centers,1)
   region = labels==k;
   region = imfill(region,'holes');
   region = 1 - region;
   labels(region==0) = k;
end
end

