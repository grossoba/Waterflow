clear all; %close all;

str_Coord = 'raster.asc';
M_Coord=importdata(str_Coord,' ',7);
grid = M_Coord.data;
xllcorner = 291564.294;
yllcorner = 63499.633;
cellsize = 1;
NODATA_value = -9999;

s_x = length(grid(1,:));
s_y = length(grid(:,1));

str_flowpath = 'flowpath.csv';
flow_Coord = importdata(str_flowpath,' ');
flowpath_x_ = flow_Coord(1:2:end,1)-xllcorner-0.5;
flowpath_y_ = flow_Coord(2:2:end,1)-yllcorner+0.5;



%% connect the points of the flow path
incr = 1;
for i = 1:length(flowpath_x_)-1
    Q1 = [flowpath_x_(i) s_y-flowpath_y_(i)];
    Q2 = [flowpath_x_(i+1) s_y-flowpath_y_(i+1)];
    %invert the points in order to manage every relative positions of the points
    if(flowpath_x_(i) > flowpath_x_(i+1))
        A = flowpath_x_(i+1);
        B = flowpath_x_(i);
    else
        A = flowpath_x_(i);
        B = flowpath_x_(i+1);
    end
    if((s_y-flowpath_y_(i)) > (s_y-flowpath_y_(i+1)))
        C = (s_y-flowpath_y_(i+1));
        D = (s_y-flowpath_y_(i));
    else
        C = (s_y-flowpath_y_(i));
        D = (s_y-flowpath_y_(i+1));
    end
    
    for x = A:B
        for y = C:D
            P = [x y];
            d = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
            
            if (d<cellsize)
                flowpath_x(incr,1) = x;
                flowpath_y(incr,1) = s_y-y;
                incr = incr + 1;
            end
        end
    end
end

% figure;
% plot(flowpath_x, s_y-flowpath_y,'xr')
% hold on;
% plot(flowpath_x_,s_y-flowpath_y_,'ob');
water = zeros(s_y,s_x);


%% erase the -9999 values
incr = 1;
for i = 1:s_x
    for j = 1:s_y
        if(grid(j,i) == -9999)
            grid(j,i) = NaN; 
        end
    end
end

%% build the initial flow path
threshold = -0.02; % -0.02 looks like a good value 
for i = 1:length(flowpath_x)
water(s_y-flowpath_y(i),flowpath_x(i)) = 1;
end


% figure;
% plot(grid(1,:));

%% propagate the water from the initial flow path
sum_init = nnz(water);
sum_res = nnz(water)+ 1;
v_choice = [-1 0 1];

figure;
surf(water);
while(sum_res > sum_init)
    sum_init = nnz(water);
    for i = 2:s_y-1
        for j = 2:s_x-1
            if(water(i,j) > 0)
            h_1 = grid(i,j)-grid(i+1,j-1);
            h_2 = grid(i,j)-grid(i+1,j); 
            h_3 = grid(i,j)-grid(i+1,j+1);
                        
            h_neigh = [h_1 h_2 h_3];
%             v_threshold = [(1-(h_1 > threshold))*(-1000)+1 (1-(h_2 > threshold))*(-1000)+1 (1-(h_3 > threshold))*(-1000)+1];
%             v_threshold = [(1-((h_1 > threshold)*(water(i+1,j-1)==0)))*(-1000)+1 (1-((h_2 > threshold)*(water(i+1,j)==0)))*(-1000)+1 (1-((h_3 > threshold)*(water(i+1,j+1)==0)))*(-1000)+1];
% 
%             if(sum(v_threshold) ~= -2997)
%             [~,index_h] = max(abs(v_threshold).*h_neigh);
%             water(i+1,j+v_choice(1,index_h)) = 2; 
%             end
            v_threshold = [];
            
            for k = 1:length(h_neigh)
                
                if((h_neigh(k) > threshold) && (water(i+1,j+v_choice(1,k)) ==0 ))
                     v_threshold(1,k) = h_neigh(k);
                else
                     v_threshold(1,k) = -999 ; 
                end
            end
            
            % test display
            if(sum(v_threshold) ~= -2997)
            h_neigh
            [water(i+1,j+v_choice(1,1)) water(i+1,j+v_choice(1,2)) water(i+1,j+v_choice(1,3))]
            v_threshold
            end 
            %
            
               if(sum(v_threshold) ~= -2997)
               index_h = 0;
               [~,index_h] = max(v_threshold);
               index_h
               water(i+1,j+v_choice(1,index_h)) = 2; 
                
               display('-----');
               end


            end
        end
    end
    sum_res = nnz(water);
end

% for i = 1:s_y
%     for j = 1:s_x
%         if(isnan(grid(i,j)))
%             water(i,j) = -1;
%         end
%     end
% end

incr_a = 1;
incr_b = 1;
for i = 1:s_y
    for j = 1:s_x
        if(isnan(grid(i,j)) || water(i,j) == 0)
            coord_water_a(incr_a,1:3) = [j i NaN];
            coord_water_b(incr_b,1:3) = [j i NaN];
            incr_a = incr_a + 1;
            incr_b = incr_b +1;
        end
        if(water(i,j) == 1)
            coord_water_a(incr_a,1:3) = [j i grid(i,j)];
            incr_a = incr_a + 1;
        end
        if(water(i,j) == 2)
            coord_water_b(incr_b,1:3) = [j i grid(i,j)];
            incr_b = incr_b + 1;
        end
    end 
end

figure;
surf(grid)
colormap summer;
hold on;
plot3(coord_water_a(:,1),coord_water_a(:,2),coord_water_a(:,3),'xb')
plot3(coord_water_b(:,1),coord_water_b(:,2),coord_water_b(:,3),'xr')

% figure;
% meshz(water)
% 
% figure;
% surf(grid)
% colormap summer;

