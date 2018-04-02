function [boundary] = frontiere( imBin )
%FRONTIERE
%
% Entrees :
%   imBin : une image binaire
%
% Sorties :
%   boundary : matrice de la frontiere de la plus grande region de imBin
%

% Calcul les frontieres de tous les objets dans imBin
[B,~] = bwboundaries(imBin);
% On trie les regions par taille
[~,ind] = sort(cellfun('length',B),'descend');
% La plus grande correspond a la forme recherchee
boundary = B{ind(1)};

end

