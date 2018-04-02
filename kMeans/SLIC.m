function [ Centers,labels ] = SLIC( im,K )
%SLIC Summary of this function goes here
%   Detailed explanation goes here

% DECOMMENTER POUR AFFICHER
% figure('Name','Image Originale')
% imshow(im);
% hold on


%% Initialisation
[nbl,nbc,c] = size(im);
imlab = rgb2lab(im);
N = nbl*nbc;
S = round(sqrt(N/K));

%% Calcul des germes
% Indices des centres
ic = round(0.5*S):S:nbl;
jc = round(0.5*S):S:nbc;
allic = repmat(ic,[length(jc) 1]);
alljc = repmat(jc,[1 length(ic)]);
% Concatenations des tous les centres initiaux
germes = [ allic(:) alljc(:) ];
%% Affinement des centres
[Glx,Gly] = gradient(imlab(:,:,1));
for k = 1:size(germes,1)
    % Calcul de voisinage 3x3 autour du germes(k)
    ivois = max(1,germes(k,1)-1):min(nbl,germes(k,1)+1);
    jvois = max(1,germes(k,2)-1):min(nbc,germes(k,2)+1);
    allivois = repmat(ivois,[length(jvois),1]);
    alljvois = repmat(jvois,[1,length(ivois)]);
    voisinage = [ allivois(:) alljvois(:) ];
    voisInd = sub2ind([nbl nbc],voisinage(:,1),voisinage(:,2));
    % Extraction des gradient du voisinage
    gradvois = Glx(voisInd) + Gly(voisInd);
    % Calcul du min
    [~,imin] = min(gradvois);
    % Remplacement du centre
    germes(k,:) = voisinage(imin,:);
end
labels = -1*ones(nbl,nbc);
distances = Inf*ones(nbl,nbc);
germAbsol = sub2ind([nbl nbc],germes(:,1),germes(:,2));
imlab1 = imlab(:,:,1);
imlab2 = imlab(:,:,2);
imlab3 = imlab(:,:,3);
Centers = [ imlab1(germAbsol) imlab2(germAbsol) imlab3(germAbsol) germes ];

%% Boucle
seuil = 10;
m = 10;
E = seuil+10;
while E > seuil
    for k = 1:size(Centers,1)
        Ck = Centers(k,:);
        %fprintf('Recherche des pixels appartenants au cluster %d\n',k); 
        for i = max(1,Ck(4)-S):min(nbl,Ck(4)+S)
            for j = max(1,Ck(5)-S):min(nbc,Ck(5)+S)
                % Distance en lab
                dc = sqrt((imlab(i,j,1) - Ck(1)).^2 + (imlab(i,j,2) - Ck(2)).^2 + (imlab(i,j,3) - Ck(3)).^2);
                % Distance physique
                ds = sqrt((i-Ck(4)).^2 + (j-Ck(5)).^2);
                %
                D = sqrt(dc^2 + ds^2*(m/S)^2);
                if D < distances(i,j)
                    % MAJ
                    distances(i,j) = D;
                    labels(i,j) = k;
                end
            end
        end
    end
    %% Renforcement de la connexite
    %labels = renforcementConnexiteOLD(labels,0.5,Centers,N);
    labels = renforcementConnexite(labels,Centers);
    
    % Calcul des nouveaux centres de clusters
    [ Centersk ] = kmeansGregAxel( Centers, imlab, labels );
    %% Calcul de l'erreur
    E = sum(sum(sqrt((Centersk(:,4:5)-Centers(:,4:5)).^2)));
    Centers = Centersk;
    
    %% Affichage des clusters par bords
    % DECOMMENTER POUR AFFICHER
%     bords = zeros(size(im(:,:,1)));
%     for k = 1:size(Centers,1)
%         bords = min(bords + edge(labels==k),1);
%     end
%     bordsInv = 1 - bords;
%     R   = 255;  % Value in range [0, 255]
%     G   = 0;
%     B   = 0;
%     RGBbords = cat(3, bords * R, bords * G, bords * B);
%     imWithoutBorders = im.*uint8(repmat(bordsInv,[1 1 3]));
%     imWithRedBorders = imWithoutBorders+uint8(RGBbords);
%     imshow(imWithRedBorders);
%     
%     plot(Centers(:,5),Centers(:,4),'r*');
%     hold on
%     drawnow;
%     pause(0.1);
    

end

% figure('Name','Appartenance aux superpixels')
% imagesc(255*labels/max(max(labels)));
% figure('Name','Image reconstruite par rapport aux superpixels');
% imreconL = zeros(nbl,nbc);
% imreconA = imreconL;
% imreconB = imreconA;
% for k = 1:size(Centers,1)
%     imreconL(find(labels==k)) = Centers(k,1);
%     imreconA(find(labels==k)) = Centers(k,2);
%     imreconB(find(labels==k)) = Centers(k,3);
% end
% imrecon(:,:,1) = imreconL;
% imrecon(:,:,2) = imreconA;
% imrecon(:,:,3) = imreconB;
% imRGB = lab2rgb(imrecon);
% subplot(2,2,1);
% imshow(imRGB(:,:,1));
% subplot(2,2,2);
% imshow(imRGB(:,:,2));
% subplot(2,2,3);
% imshow(imRGB(:,:,3));
% subplot(2,2,4);
% imshow(imRGB);
    
end

