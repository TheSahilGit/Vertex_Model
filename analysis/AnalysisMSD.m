clear; clc; close all; 

nrun = 2; 

[~,~,v_in,inn_in,num_in, forces_in] = LoadData(1000, nrun);

ct = 1;
for it = 1000:1000:400000

    it

   [Lx, Ly, v,inn,num, forces] = LoadData(it, nrun);
   MSD(ct, :) =  Calculate_MSD(Lx,Ly,v_in, inn_in, num_in, v, inn, num);
   time(ct) = it*1e-3;
   ct = ct + 1;

end

mean_MSD = mean(MSD,2);

%%

figure("Position",[100 100 800 800])

loglog(time, mean_MSD,'o', "LineWidth",3, 'MarkerSize',20);
hold on; 
loglog(time, 5e-4*time.^2, "LineWidth",3, 'DisplayName',"t^2");
hold on; 
loglog(time, 2e-3*time, "LineWidth",3, "DisplayName",'t')


xlabel("Time");
ylabel("MSD")
set(gca, "FontSize", 32)
axis square
legend()



%% Comparison with T2


X = log(time(1:end));
Y = log(mean_MSD(1:end));
p = polyfit(X, Y, 1)


figure()
loglog(time, mean_MSD, 'o', 'MarkerSize',20);
hold on
loglog(time, exp(p(2))*time.^p(1))
axis square
legend("FontSize",18)
xlabel("time")
ylabel("MSD")
axis square
set(gca, 'FontSize', 28,'LineWidth',1);
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2);

%%



fname_T2 = sprintf('../data/T2_count.dat');
fid = fopen(fname_T2);
dum4 = fread(fid,1,'float32');
T2_count = fread(fid,100000000,'float64');
cumsum_T2 = cumsum(T2_count);

figure()
plot(cumsum_T2(20000:10000:end), DisplayName="T2", LineWidth=3)
hold on
plot(mean_MSD, DisplayName="MSD", LineWidth=3)

ax = gca;

% --- Set xticks so that each tick gets a gridline ---
ax.XTick = 0:5:2000;     % <-- choose your own spacing here
                           % (whatever spacing you want gridlines at)

grid on                   % major grid = gridlines at major ticks
ax.GridLineWidth = 1.5;

legend(["cumsum T2" "MSD"], Location="northwest")
set(gca, 'FontSize', 40)




%%

%writematrix([time' mean_MSD], 'msd.dat')



%%


% figure()
% loglog(time, mean_MSD, '-o');
% hold on
% %loglog(time, msd_nrun2, '-o')
% axis square
% %legend("FontSize",18)
% xlabel("time")
% ylabel("MSD")
% axis square
% set(gca, 'FontSize', 28,'LineWidth',1);
% set(findall(gca, 'Type', 'Line'), 'LineWidth', 2);


%%

%writematrix([time' mean_MSD msd_nrun2], 'msd.dat')


%%

function [MSD] = Calculate_MSD(Lx,Ly,v_in, inn_in, num_in, v, inn, num)
    % Calculate_MSD: Computes the Mean Squared Displacement (MSD) for particles
    %
    % Inputs:
    %   v_in  - Initial positions of particles [Nx2] (x, y)
    %   inn_in - Initial neighbors indices [NxM] (up to M neighbors per particle)
    %   num_in - Initial number of neighbors per particle [Nx1]
    %   v     - Current positions of particles [Nx2] (x, y)
    %   inn   - Current neighbors indices [NxM]
    %   num   - Current number of neighbors per particle [Nx1]
    %
    % Output:
    %   MSD values are displayed for each particle.

    % Initialize parameters
    N = Lx*Ly; % Total number of particles
    MSD = zeros(N, 1); % To store MSD for each particle

    [inside1, inside2, Boundary] = Mesh_Info(Lx,Ly);
    % Loop over each particle
    % for ii = 1:length(inside2)
    %     i = inside2(ii);
    Nc = find(num ~= 0, 1, 'last')
    for i = 1:Nc
        % Get the current neighbors' positions
        vx = v(inn(i, 1:num(i)), 1); % x-coordinates of neighbors
        vy = v(inn(i, 1:num(i)), 2); % y-coordinates of neighbors

        % Compute the center of mass of neighbors
        Cx = mean(vx);
        Cy = mean(vy);

        % Get initial neighbors' positions
        vx_in = v_in(inn_in(i, 1:num_in(i)), 1); % x-coordinates of neighbors (initial)
        vy_in = v_in(inn_in(i, 1:num_in(i)), 2); % y-coordinates of neighbors (initial)

        % Compute the initial center of mass of neighbors
        Cx_in = mean(vx_in);
        Cy_in = mean(vy_in);

        % Calculate squared displacement of the center of mass
        dx = Cx - Cx_in;
        dy = Cy - Cy_in;

        % Mean Squared Displacement
        MSD(i) = dx^2 + dy^2;
    end

    % Display the MSD for each particle
    % disp('Mean Squared Displacement for each particle:');
    % disp(MSD);
end


function [inside1, inside2,Boundary] = Mesh_Info(Lx,Ly)
mainarea=(1:Lx*Ly);
leftpanel=(1:Ly);
rightpanel=(Lx*Ly-Ly+1:Lx*Ly);
toppanel=(Ly:Ly:Lx*Ly);
bottompanel=(1:Ly:Lx*Ly-Ly+1);
corners=[1 Ly Lx*Ly-Ly+1 Lx*Ly ];

leftpanel2=(Ly+1:2*Ly);
rightpanel2=(Lx*Ly-2*Ly+1:Lx*Ly-2*Ly+Ly);
toppanel2=(Ly-1:Ly:Lx*Ly-1);
bottompanel2=(2:Ly:Lx*Ly-Ly+2);
%corners2=[1 Ly Lx*Ly-Ly+1 Lx*Ly ];


leftpanel3=(2*Ly+1:3*Ly);
rightpanel3=(Lx*Ly-3*Ly+1:Lx*Ly-3*Ly+Ly);
toppanel3=(Ly-2:Ly:Lx*Ly-2);
bottompanel3=(3:Ly:Lx*Ly-Ly+3);
% corners2=[1 L L*L-L+1 L*L ];

Boundary=[leftpanel rightpanel toppanel bottompanel];

Boundary2 =[Boundary leftpanel2 rightpanel2 toppanel2 bottompanel2];

Boundary3 =[Boundary2 leftpanel3 rightpanel3 toppanel3 bottompanel3];

inside1 = setdiff(mainarea,Boundary);
inside2 = setdiff(mainarea,Boundary2);


end

