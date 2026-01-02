clear; clc; close all; 


nrun = 2; 

ct = 1;

for it = 1000:1000:2000000

[~, ~, v, inn, num, ~, biochemdata] = LoadData(it, nrun);
Nc = find(num ~= 0, 1, 'last');
%Nc = 1 ; 
Rho = biochemdata(:,1);
ROCK = biochemdata(:,2);
Myosin = biochemdata(:,3);

Rho_time(ct) = mean(Rho(1:Nc));
ROCK_time(ct) = mean(ROCK(1:Nc));
Myosin_time(ct) = mean(Myosin(1:Nc));
Myosin_time_std(ct) = std(Myosin(1:Nc));

Rho_time1(ct) = Rho(1);
ROCK_time1(ct) = ROCK(1);
Myosin_time1(ct) = Myosin(1);

[areaE,~] = AverageAreaPerimeter(Nc,v,inn,num);
%area(ct) =  areaE(1);
area(ct) = mean(areaE);
time(ct) = it * 2.5d-4;

area_individual(1, ct) = areaE(randi(Nc));
area_individual(2, ct) = areaE(randi(Nc));
area_individual(3, ct) = areaE(randi(Nc));
area_individual(4, ct) = areaE(randi(Nc));
area_individual(5, ct) = areaE(randi(Nc));




ct = ct + 1

end

%writematrix([time' area'], "mean_area_time.dat")



%%  Plots

figure('Position',[100 100 1600 800])
%figure()

subplot(1,3,1)
plot(Rho_time, 'LineWidth', 4, 'DisplayName', 'Rho')
hold on
plot(ROCK_time,   'LineWidth', 4, 'DisplayName', 'ROCK');
hold on
plot(Myosin_time, '-o', 'LineWidth', 4, 'DisplayName', 'Myosin', 'Color',[0.9 0.69 0.12])
hold on

% interval = 10;                  % plot error bar every 20 points
% idx = 1:interval:length(time);  % indices where error bars appear
% errorbar(idx, Myosin_time(idx), Myosin_time_std(idx), ...
%     'o', ...
%     'LineWidth', 2, ...      % thinner error bar lines
%     'MarkerSize', 6, ...     % smaller markers
%     'DisplayName', 'Myosin', 'Color',[0.9 0.69 0.12])



xlabel('Iteration')
ylabel('Mean Concentration')
%subtitle("No coupling")
set(gca, 'FontSize', 20, 'FontName', 'Serif')
legend('Location','northeast', 'FontSize', 20)
axis square

subplot(1,3,2)
plot(area, 'LineWidth', 4,'LineStyle', '--', 'DisplayName', 'Area')
yline(1)


xlabel('Iteration')
ylabel('Mean Area')
%subtitle("No coupling")
set(gca, 'FontSize', 20, 'FontName', 'Serif')
legend('Location','northeast', 'FontSize', 20)
axis square


subplot(1,3,3)
plot(Myosin_time,  Rho_time , 'LineWidth', 4, 'DisplayName', 'Rho')
hold on
scatter(Myosin_time(1), Rho_time(1),100,  ...
    'MarkerFaceColor','r','MarkerEdgeColor','r', 'DisplayName','Start')
hold on
scatter(Myosin_time(end), Rho_time(end),100,  ...
    'MarkerFaceColor','g','MarkerEdgeColor','g', 'DisplayName','End')
axis square



xlabel('Mean Myosin')
ylabel('Mean Rho')
%subtitle("No coupling")
set(gca, 'FontSize', 20, 'FontName', 'Serif')
legend('Location','northeast', 'FontSize', 20)


figure;

plot(area_individual(1,:), 'LineWidth', 4, 'DisplayName', 'cell 1');
hold on;
plot(area_individual(2,:), 'LineWidth', 4, 'DisplayName', 'cell 2');
plot(area_individual(3,:), 'LineWidth', 4, 'DisplayName', 'cell 3');   
hold on;
plot(area_individual(4,:), 'LineWidth', 4, 'DisplayName', 'cell 4');
plot(area_individual(5,:), 'LineWidth', 4, 'DisplayName', 'cell 5');

xlabel('time')
ylabel('area')
legend('Location','best')
set(gca, 'FontSize', 32)
grid on;



%%


function [area,perimeter] = AverageAreaPerimeter(Nc,v,inn,num)

for i = 1:Nc
    % Get the vertices for the current polygon
    vx = v(inn(i, 1:num(i)), 1);
    vy = v(inn(i, 1:num(i)), 2);

    % Calculate area
    areaE(i) = polyarea(vx, vy);

    perimeterE(i) = calculatePerimeter(vx, vy);

end
area = areaE; perimeter = perimeterE;
%area = mean(areaE); perimeter = mean(perimeterE);
end

function perimeter = calculatePerimeter(vertices_x, vertices_y)
% Ensure that the number of x and y vertices match
assert(length(vertices_x) == length(vertices_y), 'Number of x and y vertices must be the same');

% Number of vertices
num_vertices = length(vertices_x);

% Initialize perimeter
perimeter = 0;

% Calculate distance between consecutive vertices
for i = 1:num_vertices
    % Calculate index of next vertex considering cyclic boundary
    next_index = mod(i, num_vertices) + 1;

    % Calculate distance between current and next vertices
    dx = vertices_x(next_index) - vertices_x(i);
    dy = vertices_y(next_index) - vertices_y(i);

    % Add distance to perimeter
    perimeter = perimeter + sqrt(dx^2 + dy^2);
end
end








