function [ masks ] = seuillage( im, nb_images, labels, centers )
%SEGMENTATION BINAIRE
% Fonction permettant de realiser le seuillage du dinosaure
%
% Entrees :
%   im : l'image courante à binariser.
%   Centers : structure de données contenant les centres des superpixels de l'image
%	      courante.
%   labels : matrice de la taille de l'image contenant l'appartenance de
%	     chaque pixel à un superpixel donné.
%
% Sorties :
%   masks : image seuillé ne contenant que des 0 et des 1.
%

[nb_lignes, nb_colonnes, ~, ~] = size(im);
l_seuil = zeros(nb_images,1);

coefficient = 20;   % Coefficient rajouté pour améliorer le seuillage

% Calcule les coefficients de seuillage de toutes les images.
for i = 1:nb_images
    l_seuil(i) = mean(centers(:,3,i))+coefficient;
end

% Initialisation de la structure contenant les masques finaux
masks = zeros(nb_lignes, nb_colonnes, nb_images);

for ind_im = 1:nb_images
    im_seuil = ones(nb_lignes, nb_colonnes);
    
    % Realise le seuillage dans une image intermediaire
    for i = 1:size(centers,1)
        im_seuil(find(labels(:,:,ind_im) == i)) = centers(i,3,ind_im) > l_seuil(ind_im);
    end
    % Affectation du seuillage à la structure finale
    masks(:,:,ind_im) = im_seuil;
    
    % A DECOMMENTER POUR AFFICHER
    % Permet d'afficher une figure qui contient un resumé des canaux LAB et RGB du
    % dinosaure courant
    % affichage_seuil;
end

save masks;

end

