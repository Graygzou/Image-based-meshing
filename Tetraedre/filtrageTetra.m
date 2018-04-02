function [ output_args ] = filtrageTetra( triDelaunay )
% pour chaque tetraedre issu de la triangulation de Delaunay :
% Etape (1) :
% Calculer le barycentre de ses 4 sommets (on pourra utiliser plusieurs barycentres avec differents
% poids de maniere a ce que le barycentre se rapproche successivement de chaque sommet) ;


% Etape (2) : 
% Projeter chaque barycentre dans les 36 masques ;



% (3) Si un barycentre ne se projette pas dans la region du dinosaure dans au moins une des 36 images
% (im mask vaut 1), retirer le tetraedre. Attention ! il est necessaire de verifier que le barycentre
% est bien projete dans l’image avant de tester si im mask vaut 1.

end

