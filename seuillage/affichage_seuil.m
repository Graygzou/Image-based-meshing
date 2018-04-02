% --------------------------------------------------------------
% Scripts permettant d'afficher les canaux RGB, canaux LAB et
% d'autres informations concernant une image.
% --------------------------------------------------------------

% Affichage des clusters par bords
bords = zeros(nb_lignes, nb_colonnes);
for k = 1:size(centers,1)
    bords = bords + edge(labels(:,:,ind_im)==k);
end
imreconL = zeros(nb_lignes, nb_colonnes);
imreconA = imreconL;
imreconB = imreconA;
for k = 1:size(centers,1)
    imreconL(find(labels(:,:,ind_im)==k)) = centers(k,1,ind_im);
    imreconA(find(labels(:,:,ind_im)==k)) = centers(k,2,ind_im);
    imreconB(find(labels(:,:,ind_im)==k)) = centers(k,3,ind_im);
end

% Declaration de la figure
figure('Name','Ensemble d''informations utiles')

%% 1er ligne de la figure = canaux LAB
subplot(3,3,1);
imagesc(imreconL); axis equal; axis off;
subplot(3,3,2); 
imagesc(imreconA); axis equal; axis off;
subplot(3,3,3);
imagesc(imreconB); axis equal; axis off;

imrecon(:,:,1) = imreconL;
imrecon(:,:,2) = imreconA;
imrecon(:,:,3) = imreconB;
imRGB = lab2rgb(imrecon);

%% 2eme ligne de la figure = canaux RGB
subplot(3,3,4)
imshow(imRGB(:,:,1));
subplot(3,3,5)
imshow(imRGB(:,:,2));
subplot(3,3,6)
imshow(imRGB(:,:,3));

%% 3eme ligne de la figure = Informations utiles
subplot(3,3,7)
imshow(imRGB);
subplot(3,3,8)
imshow(bords);
hold on
plot(Centers(:,5,ind_im),Centers(:,4,ind_im),'r*');
subplot(3,3,9)
imshow(masks(:,:,ind_im));
