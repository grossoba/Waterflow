clear all; 
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The idea of this code is to deform a hexagonal grid adding points that
% will attract or repulse the vertices. We try to keep the hexagonal
% structure but for extreme cases it can break.
%
% The parts of the code are organized this way : 
%   1) define the spacing distance and the size of the system
%   2) create the hexagonal grid structure
%   3) define the influence points (as many as desired) with in this order 
%      (coordinate x,coordinate y,sigma,force); the repulsion or attraction
%      is given by the sign of "force" and sigma defines the inluence zone
%      around the influence point. The transformation of the coordinates is
%      then done this way : coord' = coord + f*exp(-dist^2/2sigma), where
%      dist is the distance between the influence point and the vertex
%      considered. In order to get a smoother deformation one can tune it 
%      with sigma and force. 
%   4) deformation of the grid according to the structures defined above 
%   5) plot of the result (first figure is the deformation field and
%      second figure is the final result)
%
%  (Bastien Grosso, February 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1)Initial Geometry (Set the initial geometry)
% (CHANGE TO MODIFY THE SIZE OF THE GRID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 4; % spacing parameter
N_hex = 20; % vertical
Nrep = 20; % horizontal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2)Creates the initial lattice (nothing to change here)
% (NOTHING TO CHANGE HERE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Coord = zeros(N_hex,2);
str_Nrep = int2str(Nrep);
N = N_hex*4 + 2; % number of atoms
a1 = [3*a,0];
Coord(1,1) = a/2;
Coord(1,2) = 0;
Coord(2,1) = 3*a/2;
Coord(2,2) = 0;
Coord(3,1) = 0;
Coord(3,2) = sqrt(3)/2*a;
Coord(4,1) = 2*a;
Coord(4,2) = sqrt(3)/2*a;

for i = 5:N
    Coord(i,1) = Coord(i-4,1);
    Coord(i,2) = Coord(i-4,2)+a*sqrt(3);
end

for i = 1:N
    for j = 1:Nrep
        Coord(i+N*j,1) = Coord(i,1)+3*a*j;
        Coord(i+N*j,2) = Coord(i,2);
    end
end

middle_x = (max(Coord(:,1)) + min(Coord(:,1)))/2;
middle_y = (max(Coord(:,2)) + min(Coord(:,2)))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3)Influence points (define here the points with [x_coord, y_coord, z_coord
% or a value, the influence (sign and intensity)])
% (CHANGE TO ADD THE INFLUENCE POINTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Influence_points = zeros(1,4);
Influence_points(1,:) = [middle_x/3,middle_y/3,80,-5];
Influence_points(2,:) = [6*middle_x/4,6*middle_y/4,50,8];
Influence_points(3,:) = [middle_x/3,6*middle_y/4,50,8];
Influence_points(4,:) = [6*middle_x/4,middle_y/3,20,4];
Influence_points(5,:) = [middle_x-2.5,middle_y+1,60,-5];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4)Dynamics (moves the grid according to the influence points) 
% and set the neighbors
% (NOTHING TO CHANGE HERE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutoff = 4000;
Displ = zeros(length(Coord(:,1)),2);
for i = 1:length(Coord(:,1))
    for j = 1:length(Influence_points(:,1))
        vect(i,:) = (Coord(i,:)-Influence_points(j,1:2));
        dist2(i) = norm(vect(i,:));
        if((dist2(i) < cutoff) && (dist2(i) > 0.1))
            Displ(i,:) = Displ(i,:) + (vect(i,:)/dist2(i))*sign(Influence_points(j,4))*min(abs(Influence_points(j,4)*exp(-dist2(i)^2/(2*(Influence_points(j,3))^2))),dist2(i)*0.8);
        else
            Displ(i,:) = Displ(i,:) + [0 0];
        end
    end
    Coord_d(i,:) = Coord(i,:) + Displ(i,:);
    
    % To garantee that vertices don't go out of the boundaries
    Coord_d(i,1) = (Coord_d(i,1) <= middle_x)*max(min(Coord(:,1)),Coord_d(i,1)) + (Coord_d(i,1) > middle_x)*min(max(Coord(:,1)),Coord_d(i,1));
    Coord_d(i,2) = (Coord_d(i,2) <= middle_y)*max(min(Coord(:,2)),Coord_d(i,2)) + (Coord_d(i,2) > middle_y)*min(max(Coord(:,2)),Coord_d(i,2));
    % To fix the boundaries
    Coord_d(i,1) = (Coord(i,1) == max(Coord(:,1)) || Coord(i,1) == min(Coord(:,1)))*Coord(i,1) + (Coord(i,1) ~= max(Coord(:,1)) && Coord(i,1) ~= min(Coord(:,1)))*Coord_d(i,1);
    Coord_d(i,2) = (Coord(i,2) == max(Coord(:,2)) || Coord(i,2) == min(Coord(:,2)))*Coord(i,2) + (Coord(i,2) ~= max(Coord(:,2)) && Coord(i,2) ~= min(Coord(:,2)))*Coord_d(i,2);

end

% Neighbors
bond1 = ones(size(Coord,1),3)*(0); % matrix that tells the first neighbor for each atom
nbond1 = zeros(1,size(Coord,1)); % number of first neighbors

for j = 1:size(Coord,1)
    for m = 1:size(Coord,1)
        
        dist2 = sum((Coord(j,:)-Coord(m,:)).^2);
        if(j~=m && dist2<(a+0.2)^2) % first neighbor
            nbond1(j) = nbond1(j) + 1;
            bond1(j,nbond1(j)) = m;
        end
    end
end

%%%%%%%%%%%%
% 5)Plots
%%%%%%%%%%%%
figure;
hold on ;
plot(Coord(:,1),Coord(:,2),'xg')
plot(Coord_d(:,1),Coord_d(:,2),'ob')
plot(Influence_points(:,1),Influence_points(:,2),'xr');
quiver(Coord(:,1),Coord(:,2),Displ(:,1),Displ(:,2))
axis equal;
axis off;

figure;
hold on;
axis off;
plot(Influence_points(:,1),Influence_points(:,2),'xr');
for i = 1:length(bond1(:,1))
    for j = 1:nbond1(i)
    plot([Coord_d(i,1) Coord_d(bond1(i,j),1)],[Coord_d(i,2) Coord_d(bond1(i,j),2)],'-b');
    end
end

    
    
    
    
    
    
    
    