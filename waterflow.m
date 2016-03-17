clear all; close all;

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

figure;
plot(flowpath_x, s_y-flowpath_y,'xr')
hold on;
plot(flowpath_x_,s_y-flowpath_y_,'ob');
water = zeros(s_y,s_x);
water_init = water;

%% erase the -9999 values
incr = 1;
for i = 1:s_x
    for j = 1:s_y
        
        if(grid(j,i) == -9999)
            Coord(incr,1:3) = [cellsize*i cellsize*j NaN];
            grid(j,i) = NaN;
        else
            Coord(incr,1:3) = [cellsize*i cellsize*j grid(j,i)];
        end
        
        incr = incr + 1;
    end
end

%% build the flow path
threshold = -0.02; % -0.02 looks like a good value 
neighbors = 1;
for i = 1:length(flowpath_x)
water(s_y-flowpath_y(i),flowpath_x(i)) = 1;
end

water_init = water;

sum_init = sum(sum(water));
sum_res = sum(sum(water))+ 1;
while(sum_res > sum_init)
    sum_init = sum(sum(water));
    for i = 2:s_y-1
        for j = 2:s_x-1
            
            if(water(i,j) == 1 && neighbors == 1)
                if(grid(i,j)-grid(i+1,j) > threshold)
                    water(i+1,j) = 1;
                end
                if(grid(i,j)-grid(i-1,j) > threshold)
                    water(i-1,j) = 1;
                end
                if(grid(i,j)-grid(i,j+1) > threshold)
                    water(i,j+1) = 1;
                end
                if (grid(i,j)-grid(i,j-1) > threshold)
                    water(i,j-1) = 1;
                end
            end
            if(water(i,j) == 1 && neighbors == 2)
                if(grid(i,j)-grid(i+1,j+1) > threshold)
                    water(i+1,j+1) = 1;
                end
                if(grid(i,j)-grid(i+1,j-1) > threshold)
                    water(i+1,j-1) = 1;
                end
                if(grid(i,j)-grid(i-1,j+1) > threshold)
                    water(i-1,j+1) = 1;
                end
                if (grid(i,j)-grid(i-1,j-1) > threshold)
                    water(i-1,j-1) = 1;
                end
            end
        end
    end
    sum_res = sum(sum(water));
end

for i = 1:s_y
    for j = 1:s_x
        if(isnan(grid(i,j)))
            water(i,j) = 2;
        end
    end
end

incr = 1;
for i = 1:s_y
    for j = 1:s_x
        if(isnan(grid(i,j)) || water(i,j) == 0)
            coord_water(incr,1:3) = [j i NaN];
            incr = incr + 1;
        end
        if(water(i,j) == 1)
            coord_water(incr,1:3) = [j i grid(i,j)];
            incr = incr + 1;
        end
    end 
end

figure;
surf(grid)
colormap summer;
hold on;
plot3(coord_water(:,1),coord_water(:,2),coord_water(:,3),'xb')

figure;
meshz(water)

