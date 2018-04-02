clear
close all

projectpath = genpath(pwd);
addpath(projectpath);

load donnees;
%Calcul des faces du maillage a garder
FACES = [];
%% Construction de toutes les faces
for i = 1:size(tetrasave,1)
    FACES = [ FACES ;
              sort([tetrasave(i,1) tetrasave(i,2) tetrasave(i,3)]) ; ...
              sort([tetrasave(i,2) tetrasave(i,3) tetrasave(i,4)]) ; ...
              sort([tetrasave(i,3) tetrasave(i,4) tetrasave(i,1)]) ; ...
              sort([tetrasave(i,4) tetrasave(i,1) tetrasave(i,2)]) ; ...
              ];
end
fprintf('Nombre de faces : %d\n',size(FACES,1));
FACES = sortrows(FACES);
fprintf('Suppression des faces qui appartiennent � plusieurs t�tra�dres\n');
i = 1;
while (i < size(FACES,1))
    FaceCour = FACES(i,:);
    if (FaceCour == FACES(i+1,:))
        itrouver = i;
        while (FaceCour == FACES(itrouver,:) & itrouver < size(FACES,1))
            itrouver = itrouver + 1;
        end
        FACES(i:itrouver,:) = [];
    else
        i = i + 1;
    end
end
fprintf('Calcul du maillage final termine : %d faces. \n',size(FACES,1));

%Affichage du maillage final
figure('Name','Maillage final de l''objet en 3D');
hold on
for i = 1:size(FACES,1)
   plot3([X(1,FACES(i,1)) X(1,FACES(i,2))],[X(2,FACES(i,1)) X(2,FACES(i,2))],[X(3,FACES(i,1)) X(3,FACES(i,2))],'r');
   plot3([X(1,FACES(i,1)) X(1,FACES(i,3))],[X(2,FACES(i,1)) X(2,FACES(i,3))],[X(3,FACES(i,1)) X(3,FACES(i,3))],'r');
   plot3([X(1,FACES(i,3)) X(1,FACES(i,2))],[X(2,FACES(i,3)) X(2,FACES(i,2))],[X(3,FACES(i,3)) X(3,FACES(i,2))],'r');
end;
