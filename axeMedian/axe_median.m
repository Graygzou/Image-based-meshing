function [ A, coordonnees ] = axe_median( Pts, Pts_bords, nb_pts_bords )
% Fonction permettant de réaliser l'axe median
%% Inputs
% ind_im : n'existera plus.
% Pts : Ensemble des points calculés à partie de Voronoi
% Pts_bords : Ensemble des points du bords grâce à la fonction Matlab
% nb_pts_bords : Nombre de points du bords.

%% Ouputs
% A : matrice d'adjacence
% A : matrice d'adjacence

C = [];
Pts = Pts';
% ATTENTION /!\
% On filtre les points qui se superpose
[C,~,~] = unique(Pts,'rows');
Pts = C';

nb_pts = size(Pts,2);

% Declaration de la matrice d'adjacence
A = zeros(nb_pts,nb_pts);

dist = zeros(nb_pts);
dist_bords = zeros(nb_pts,nb_pts_bords);

% Calcul les distances eucliennnes( entre tous les points
for i = 1:nb_pts
    % Recupere le point courant
    pts_courant = Pts(:,i);

    % Calcul les distances euclidiennes entre CE POINT et chaque AUTRE points
    dist(i,:) = sqrt(sum((repmat(pts_courant,1,nb_pts) - Pts).^2,1));
    
    % Calcul les distances euclidiennes entre CE point et TOUS les points du bord. 
    dist_bords(i,:) = sqrt(sum((repmat(pts_courant,1,nb_pts_bords) - Pts_bords).^2,1));
end
Rayons =  min(dist_bords,[],2);
Rayons = Rayons';

% ATTENTION /!\
% On filtre les points abérrants qui n'appartiennent pas au dinosaure
[~, ind_valide] = find(Rayons < 100);
Rayons = Rayons(ind_valide);
Pts = Pts(:,ind_valide);
dist = dist(ind_valide,ind_valide);

% Initialise les points que l'on retourna.
coordonnees = [Pts(2,:)', Pts(1,:)'];

% Verification avec l'affichage des cercles
% for u = 1:size(Pts,2)
%     t = 0:360;
%     t = t * pi/ 180;
%     ind_cercle_x = Pts(1,u) + cos(t) * Rayons(u);
%     ind_cercle_y = Pts(2,u) + sin(t) * Rayons(u);
% 
%     hold on;
%     plot(Pts(2,u), Pts(1,u),'m*')
%     plot(ind_cercle_y, ind_cercle_x)
%     hold off;
% end


%% Debut de l'algo
ind_dispo = [];              % Indices du cluster courant.
seuil = 25;           % seuil autorisé pour qu'un point soit candidat.
decroissance = 0.1;     % coefficient de decroissance du rayon.
rayon_courant = max(Rayons);  % Commence avec le cercle le plus grand.

while ~(all(Rayons > rayon_courant))
                    
    % On recherche les points dont le cercle est >= rayon_courant
    ind_trouve = find(Rayons >= rayon_courant);
    % on filtre les indices deja traité (\in ind_dispo)
    for u = ind_dispo
        ind_trouve = ind_trouve(find(ind_trouve ~= u));
    end
    
    % Verification affichage des cercles des points trouvés
%     for u = ind_trouve
%         t = 0:360;
%         t = t * pi/ 180;
%         ind_cercle_x = Pts(1,u) + cos(t) * Rayons(u);
%         ind_cercle_y = Pts(2,u) + sin(t) * Rayons(u);
% 
%         hold on;
%         h2 = plot(Pts(2,u), Pts(1,u),'b.');
%         h3 = plot(ind_cercle_y, ind_cercle_x,'Color','b');
%         hold off;
%     end
%     pause(0.001);

    if length(ind_trouve) ~= 0
        % Si c'est le premier point
        if length(ind_dispo) < 1
            
            % On le rajoute aux pts dispo
            ind_dispo = [ind_dispo, ind_trouve];
            ind_trouve = [];
        else
            ind_trouve;
            ind_dispo;

            % on filtre les indices qui n'appartiennent pas a un certain seuil
            % autour du cluster
            dist_cluster = zeros(1,length(ind_trouve));
            for i = 1:length(ind_trouve)
                dist_cluster(i) = min(dist(ind_trouve(i),ind_dispo));
            end
            ind_trouve = ind_trouve(find(dist_cluster <= seuil));

            % tant qu'il reste des points a accrocher au cluster
            while length(ind_trouve) > 0
                dist_cluster = zeros(1,length(ind_trouve));
                % Trie les points trouve par leur distance
                for i = 1:length(ind_trouve)
                    dist_cluster(i) = min(dist(ind_trouve(i),ind_dispo));
                end
                % Recupere le nouveau point plus proche du cluster
                [~, inf_dist_toCluster] = min(dist_cluster);
                pts_courant_trouve = ind_trouve(inf_dist_toCluster);

                % Cherche le point dans le cluster le plus proche du point courant
                [dist_entre2,ind_min_existant] = min(dist(pts_courant_trouve,ind_dispo));

                % On regarde s'il n'existe des points entre ceux deux
                % points
                dist_decoupe_max = sqrt((dist_entre2/2)^2 + (1.8)^2)*2;
                % trouve les points "preque" equidistant des deux points
                % bientot liés
                pts_entre2 = find((dist(pts_courant_trouve,:) + dist(ind_dispo(ind_min_existant),:)) <= dist_decoupe_max);
                
                if ~isempty(pts_entre2)
                    
                    % Filtrages les indices deja traité (ind_dispo)
                    for o = ind_dispo
                        pts_entre2 = pts_entre2(pts_entre2 ~= o);
                    end
                    pts_entre2 = pts_entre2(pts_entre2 ~= pts_courant_trouve);
                    
                    if ~isempty(pts_entre2)
                    
                        % On choisit le point le plus proche du cluster
                        [dist_outsider,ind_outsider] = min(dist(ind_dispo(ind_min_existant),pts_entre2));

                        if dist_outsider < dist_entre2
                            % Fait la liaison entre le cluster et ce nouveau point
                            A(ind_dispo(ind_min_existant), pts_entre2(ind_outsider)) = 1;
                            %A(pts_entre2(ind_outsider),ind_dispo(ind_min_existant)) = 1;
                            % Affichage
%                             hold on;
%                             line([Pts(2,pts_entre2(ind_outsider)), Pts(2,ind_dispo(ind_min_existant))], [Pts(1,pts_entre2(ind_outsider)), Pts(1,ind_dispo(ind_min_existant))], ...
%                                 'Color','red','LineWidth',2,'LineStyle','-');
%                             hold off;

                            % Mise a jour
                            ind_dispo = [ind_dispo, pts_entre2(ind_outsider)];
                        else
                            % Fait la liaison
                            A(ind_dispo(ind_min_existant), pts_courant_trouve) = 1;
                            %A(pts_courant_trouve,ind_dispo(ind_min_existant)) = 1;
                            % Affichage
%                             hold on;
%                             line([Pts(2,pts_courant_trouve), Pts(2,ind_dispo(ind_min_existant))], [Pts(1,pts_courant_trouve), Pts(1,ind_dispo(ind_min_existant))], ...
%                                 'Color','red','LineWidth',2,'LineStyle','-');
%                             hold off;

                            % Mise a jour
                            ind_dispo = [ind_dispo, pts_courant_trouve];
                            ind_trouve = ind_trouve(ind_trouve ~= pts_courant_trouve);
                        end
                    else
                        % Fait la liaison
                        A(ind_dispo(ind_min_existant), pts_courant_trouve) = 1;
                        %A(pts_courant_trouve,ind_dispo(ind_min_existant)) = 1;
                        % Affichage
%                         hold on;
%                         line([Pts(2,pts_courant_trouve), Pts(2,ind_dispo(ind_min_existant))], [Pts(1,pts_courant_trouve), Pts(1,ind_dispo(ind_min_existant))], ...
%                            'Color','red','LineWidth',2,'LineStyle','-');
%                         hold off;

                        % Mise a jour
                        ind_dispo = [ind_dispo, pts_courant_trouve];
                        ind_trouve = ind_trouve(ind_trouve ~= pts_courant_trouve);
                    end
                else
                    % Fait la liaison
                    A(ind_dispo(ind_min_existant), pts_courant_trouve) = 1;
                    %A(pts_courant_trouve,ind_dispo(ind_min_existant)) = 1;
                    %hold on;
                    %line([Pts(2,pts_courant_trouve), Pts(2,ind_dispo(ind_min_existant))], [Pts(1,pts_courant_trouve), Pts(1,ind_dispo(ind_min_existant))], ...
                    %    'Color','red','LineWidth',2,'LineStyle','-');
                    %hold off;

                    % Mise a jour
                    ind_dispo = [ind_dispo, pts_courant_trouve];
                    ind_trouve = ind_trouve(ind_trouve ~= pts_courant_trouve);
                end
            end
        end
    end

    % Mise a jour
    rayon_courant = rayon_courant - decroissance;
    
end

end

