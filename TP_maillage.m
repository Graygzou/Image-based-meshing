close all;
nb_images = 36; % Nombre d'images

projectpath = genpath(pwd);
addpath(projectpath);

%% chargement des images
for i = 1:nb_images
    if i<=10
        nom = sprintf('images/viff.00%d.ppm',i-1);
    else
        nom = sprintf('images/viff.0%d.ppm',i-1);
    end;
    % L'ensemble des images de taille : nb_lignes x nb_colonnes x nb_canaux
    % x nb_images
    im(:,:,:,i) = imread(nom); 
    
end;
% chargement des points 2D suivis
% pts de taille nb_points x (2 x nb_images)
% sur chaque ligne de pts 
% tous les appariements possibles pour un point 3D donne
% on affiche les coordonnees (xi,yi) de Pi dans les colonnes 2i-1 et 2i
% tout le reste vaut -1
pts = load('viff.xy');
% Chargement des matrices de projection
% Chaque P{i} contient la matrice de projection associee a l'image i 
% RAPPEL : P{i} est de taille 3 x 4
load dino_Ps;
% chargement des masques (pour l'elimination des fonds bleus)
% de taille nb_lignes x nb_colonnes x nb_images

[nb_lignes, nb_colonnes, nb_canaux, ~] = size(im);
if (exist('SLIC') == 2)
    load SLIC;
end

%% Segmentation en superpixels
if ~exist('Centers','var') || ~exist('labels','var')
    nb_super_pixel = input('Veuillez saisir le nombre de superpixels : ');

    % Calcul les super-pixels pour la premiere image
    % pour avoir les dimensions précises
    fprintf('Traitement de l''image 1\n');
    [Centers1, labels1] = SLIC(im(:,:,:,1), nb_super_pixel);

    % Initialise les structures
    Centers = zeros(size(Centers1,1), size(Centers1,2), nb_images);
    labels = zeros(size(labels1,1), size(labels1,2), nb_images);
    Centers(:,:,1) = Centers1;
    labels(:,:,1) = labels1;

    % Calcul les super-pixels pour les images restantes (2:end)
    for ind_im = 2:nb_images
        fprintf('Traitement de l''image %d\n',ind_im);
        [Centers(:,:,ind_im), labels(:,:,ind_im)] = SLIC( im(:,:,:,ind_im),nb_super_pixel);
    end
    
    save SLIC;
end

%% Effectue le seuillage
% permet de binariser les images
% COMMENTER LE IF SI
fprintf('----------- Debut Seuillage ---------------\n');
masks = seuillage( im, nb_images, labels, Centers );
fprintf('-------------- Fin Seuillage -----------------\n');

%% Affichage des images
figure('Name','Images'); 
subplot(2,2,1); imshow(im(:,:,:,1)); title('Image 1');
subplot(2,2,2); imshow(im(:,:,:,9)); title('Image 9');
subplot(2,2,3); imshow(im(:,:,:,17)); title('Image 17');
subplot(2,2,4); imshow(im(:,:,:,25)); title('Image 25');

%% Affichage des masques associes
figure('Name','Masques associes');
subplot(2,2,1); imshow(masks(:,:,1)); title('Masque image 1');
subplot(2,2,2); imshow(masks(:,:,9)); title('Masque image 9');
subplot(2,2,3); imshow(masks(:,:,17)); title('Masque image 17');
subplot(2,2,4); imshow(masks(:,:,25)); title('Masque image 25');

%% Reconstruction des points 3D
X = []; % Contient les coordonnees des points en 3D
color = []; % Contient la couleur associee
if nb_images == 36
% Pour chaque coupple de points apparies
for i = 1:size(pts,1)
    % Recuperation des ensembles de points apparies
    l = find(pts(i,1:2:end)~=-1);
    % Verification qu'il existe bien des points apparies dans cette image
    if size(l,2) > 1 & max(l)-min(l) > 1 & max(l)-min(l) < 36
        A = [];
        R = 0;
        G = 0;
        B = 0;
        % Pour chaque point recupere, calcul des coordonnees en 3D
        for j = l
            A = [A ; P{j}(1,:) - pts(i,(j-1)*2+1) * P{j}(3,:);
            P{j}(2,:)-pts(i,(j-1)*2+2)*P{j}(3,:)];
            R = R + double(im(int16(pts(i,(j-1)*2+1)),int16(pts(i,(j-1)*2+2)),1,j));
            G = G + double(im(int16(pts(i,(j-1)*2+1)),int16(pts(i,(j-1)*2+2)),2,j));
            B = B + double(im(int16(pts(i,(j-1)*2+1)),int16(pts(i,(j-1)*2+2)),3,j));
        end;
        [U,S,V] = svd(A);
        X = [X  V(:,end)/V(end,end)];
        color = [color  [R/size(l,2);G/size(l,2);B/size(l,2)]];
    end;
end;
fprintf('Calcul des points 3D termine : %d points trouves. \n',size(X,2));
end

%% Affichage du nuage de points 3D
figure;
hold on;
for i = 1:size(X,2)
    plot3(X(1,i),X(2,i),X(3,i),'.','col',color(:,i)/255);
end;
axis equal;

%% Estimation de l'axe median
boudaries = [];
fprintf('---------------- Debut Axe_median ------------------\n');
for ind_im = 1:nb_images
    
    % Estimation de la frontiere
    boudaries = frontiere(masks(:,:,ind_im));
    
    % Estimation des points du squelette
    [vx, vy] = voronoi(boudaries(:,1),boudaries(:,2));
    
    vx = round(vx(1,:));
    vy = round(vy(1,:));
    
    %% Filtre les points qui sont en dehors de l'objet
    indices = [];
    for i = 1 : size(vy,2)
       if vx(1,i) > 0 && vx(1,i) <= nb_lignes && vy(1,i) > 0 && vy(1,i) <= nb_colonnes
           if masks(vx(1,i), vy(1,i),ind_im) == 1
               indices = [indices i];
           end
       end
    end

    % DECOMMENTER POUR AFFICHER
    % Affichage du contour et des points issus de voronoi
%     figure('Name','Frontiere');
%     imagesc(masks(:,:,ind_im))
%     colormap gray;
%     hold on;
%     plot(boudaries(:,2),boudaries(:,1), 'g', 'LineWidth', 2)
%     plot(vy(1,indices),vx(1,indices),'b.')
%     hold off;
%     pause(0.001)
    
    %% Construction de l'axe median
    [A, coordonnees] = axe_median([vx(1,indices); vy(1,indices)], [boudaries(:,1)'; boudaries(:,2)'], length(boudaries(:,1)));
    
    % DECOMMENTER POUR AFFICHER L'AXE MEDIAN
    % Affichage l'axe median obtenu
%     figure('Name','Frontiere');
%     imagesc(masks(:,:,ind_im))
%     colormap gray;
%     hold on;
%     gplot(A, coordonnees,'r-')
%     pause(0.001);
end
fprintf('------------- Fin Axe_median -------------\n');

% Tetraedrisation de Delaunay
T = delaunayTriangulation(X(1,:)',X(2,:)',X(3,:)');

% A DECOMMENTER POUR AFFICHER LE MAILLAGE
fprintf('Tetraedrisation terminee : %d tetraedres trouves. \n',size(T,1));
% Affichage de la tetraedrisation de Delaunay
figure;
tetramesh(T);

% Calcul des barycentres de chacun des tetraedres
% Coeff qui permet de projeter en 2D ? Bizarre..
poids = [ 1/4 1/4 1/4 1/4 ; ...
          1/3 2/9 2/9 2/9 ; ...
          2/9 1/3 2/9 2/9 ; ...
          2/9 2/9 1/3 2/9 ; ...
          2/9 2/9 2/9 1/3 ];
nb_barycentres = 5;
for i = 1:size(T,1)
    % Calcul des barycentres differents en fonction des poids differents
    % En commencant par le barycentre avec poids uniformes
    indicesCour = T(i,:);
    % On calcule 5 barycentres chacun avec une pond�ration differente
    C_g(:,i,1)= [ poids(1,:)*T.Points(indicesCour,:) 1 ]';
    C_g(:,i,2)= [ poids(2,:)*T.Points(indicesCour,:) 1 ]';
    C_g(:,i,3)= [ poids(3,:)*T.Points(indicesCour,:) 1 ]';
    C_g(:,i,4)= [ poids(4,:)*T.Points(indicesCour,:) 1 ]';
    C_g(:,i,5)= [ poids(5,:)*T.Points(indicesCour,:) 1 ]';
end

% A DECOMMENTER POUR VERIFICATION 
% A RE-COMMENTER UNE FOIS LA VERIFICATION FAITE
% Visualisation pour verifier le bon calcul des barycentres
% for i = 1:nb_images
%    fprintf('Bary sur img%d\n',i);
%    for k = 1:nb_barycentres
%        fprintf('%d\n',k);
%        o = P{i}*C_g(:,:,k);
%        o = o./repmat(o(3,:),3,1);
%        imshow(im_mask(:,:,i));
%        hold on;
%        plot(o(2,:),o(1,:),'rx');
%        pause;
%        close;
%    end
%    fprintf('Fin img%d\n',i);
% end


% A DECOMMENTER ET A COMPLETER
% Copie de la triangulation pour pouvoir supprimer des tetraedres
% Retrait des tetraedres dont au moins un des barycentres 
% ne se trouvent pas dans au moins un des masques des images de travail
% A COMMENTER POUR UTILISER LES MASQUES FOURNIS
load mask.mat

masks = im_mask;
%
tetrasave = T.ConnectivityList;
C_g_bis = C_g;
% Pour chaque barycentre
for indMask = 1:size(P,2)
    % Masque courant
    maski = masks(:,:,indMask);
    for k = 1:nb_barycentres
        % Projection des barycentres
        barysK = C_g_bis(:,:,k);
        baryProj = P{indMask}*barysK;
        baryProj = round(baryProj./repmat(baryProj(3,:),3,1));
        % Verifications de la validite des coords
        if ( baryProj(1,:) > 0 & baryProj(1,:) <= size(im_mask,1) & baryProj(2,:) > 0 & baryProj(2,:) <= size(im_mask,2))
            %Indices pour le masque
            indBary = sub2ind(size(maski),baryProj(1,:),baryProj(2,:));
            % Indices des tetraedes dont les barycentres ne conviennent pas
            % au masque
            indNotOk = find(maski(indBary)==1);
            % Suppression des ces tetraedres 
            C_g_bis(:,indNotOk,:) = [];
            tetrasave(indNotOk,:) = [];
        end 
    end
end

% A DECOMMENTER POUR AFFICHER LE MAILLAGE RESULTAT
% Affichage des tetraedres restants
fprintf('Retrait des tetraedres exterieurs a la forme 3D termine : %d tetraedres restants. \n',size(tetrasave,1));
figure('Name','Maillage epure');
trisurf(tetrasave,X(1,:),X(2,:),X(3,:));

% Sauvegarde des donnees
save donnees;
