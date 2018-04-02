function [ newCenters ] = kmeansGregAxel( Centers, imlab, labels )
%KMEANSGREGAXEL 
% Calcule les nouveaux centres des superpixels en fonctions des pixels qui
% les composent.
%
% Entrees :
%   Centers : liste des centres de superpixels [ l a b x y ]
%   imlab : image considere en LAB
%   labels : labels d'appartenances de chaque pixel a un superpixel
%
% Sorties :
%   newCenters : nouveaux centres des superpixels 
%
newCenters = [];
L = imlab(:,:,1);
A = imlab(:,:,2);
B = imlab(:,:,3);
for k = 1:size(Centers,1)
    % Pour chaque centre
    % on recupere les indices des pixels appartenant au cluster k
    indCluster = find(labels==k);
    if isempty(indCluster)
        newCenters = [ newCenters ; zeros(1,5) ];
    else
        [xCluster,yCluster] = ind2sub(size(L),indCluster);
        % On recupere les donnees de ces pixels : l a b x y
        imCluster = [ L(indCluster) A(indCluster) B(indCluster) xCluster yCluster ];
        % Le nouveau centre du cluster est egal a la moyenne des pixels du
        % cluster
        newCenters = [ newCenters ; mean(imCluster,1) ];
    end
end
newCenters(:,4:5) = round(newCenters(:,4:5));

end

