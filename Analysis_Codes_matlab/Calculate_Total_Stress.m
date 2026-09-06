function [total_stress_tensor, ShearStress_Individual] = Calculate_Total_Stress(Lx,Ly,v,inn,num,radius, beta_arr)

% -------------------- PARAMETERS --------------------
% Was reading these via readtable()+magic row numbers (para1(1,1) etc.),
% which only worked by relying on an undocumented MATLAB readtable quirk
% (see ReadPara1Params.m's header) -- switched to the name-based reader so
% this stays correct if para_Simulation.dat's line count ever changes
% again (it just did, for if_PBC).
p1 = ReadPara1Params("../para_Simulation.dat");
A0     = p1.Ao;
lambda = p1.lambda;

% BUGFIX (log.txt, compute_CellStress.m): p1.beta/p1.gamm are the RAW
% para_Simulation.dat values. allocation.f90::read_input nondimensionalizes
% them exactly ONCE right after reading (beta = beta/(lambda*Ao); gamm =
% gamm/(lambda*Ao^1.5)) before Force.f90/Stress.f90 ever use them -- so the
% values the simulation's own stress tensor actually uses are the RESCALED
% ones, not the file's raw numbers. Was previously read raw here, silently
% wrong whenever lambda*Ao (or lambda*Ao^1.5) isn't 1 -- a no-op for this
% repo's current para_Simulation.dat (lambda=Ao=1) but wrong in general.
beta   = p1.beta / (lambda * A0);
gamma  = p1.gamm / (lambda * A0^1.5);

% Optional beta array handling
if nargin < 7 || isempty(beta_arr)
    use_beta_array = false;
else
    use_beta_array = true;
end

% -------------------- DOMAIN INFO --------------------
para2 = load(strcat("../para_MeshDims.dat"));
Lx = para2(1);
Ly = para2(2);

% -------------------- CELL SELECTION --------------------

[inside_cells, ~, ~, ~] = Find_Cells_Within_Radius(Lx, Ly, v, inn, num, radius);
inside2 = inside_cells; 

% -------------------- INITIALIZATION --------------------
TotalArea = 0;
total_stress_tensor = zeros(2);
ShearStress_Individual = zeros(length(inside2),1);

% -------------------- TOTAL AREA --------------------
for ii = 1:length(inside2)
    ic = inside2(ii);

    vx = v(inn(ic,1:num(ic)), 1);
    vy = v(inn(ic,1:num(ic)), 2);
    [vx, vy] = unwrapCellPBC(vx, vy, Lx, Ly);

    [geom, ~, ~] = polygeom(vx, vy);
    area = abs(geom(1));

    TotalArea = TotalArea + area;
end

% -------------------- STRESS COMPUTATION --------------------
for ii = 1:length(inside2)

    ic = inside2(ii);

    % Reset sigma per cell (IMPORTANT FIX)
    sigma = zeros(2);

    vx = v(inn(ic,1:num(ic)), 1);
    vy = v(inn(ic,1:num(ic)), 2);
    [vx, vy] = unwrapCellPBC(vx, vy, Lx, Ly);

    [geom, ~, ~] = polygeom(vx, vy);
    area = abs(geom(1));
    perimeter = abs(geom(4));

    % --- Beta selection ---
    if use_beta_array
        beta_local = beta_arr(ic);
    else
        beta_local = beta;
    end

    % --- Terms ---
    term1 = 2.0d0 * lambda * (area - A0);
    term2 = (2.0d0 * beta_local * perimeter + gamma)/(2.0d0*area);

    % --- Edge loop ---
    for j = 1:num(ic)
        jp = j + 1;
        if j == num(ic)
            jp = 1;
        end

        X = vx(jp) - vx(j);
        Y = vy(jp) - vy(j);
        R = sqrt(X*X + Y*Y);

        sigma(1,1) = sigma(1,1) + X*X/R;
        sigma(1,2) = sigma(1,2) + X*Y/R;
        sigma(2,1) = sigma(2,1) + Y*X/R;
        sigma(2,2) = sigma(2,2) + Y*Y/R;
    end

    % --- Final stress tensor for cell ---
    sigma(1,1) = term1 + term2 * sigma(1,1);
    sigma(1,2) = term2 * sigma(1,2);
    sigma(2,1) = term2 * sigma(2,1);
    sigma(2,2) = term1 + term2 * sigma(2,2);

    % --- Store shear contribution ---
    ShearStress_Individual(ii) = sigma(1,2) * area / TotalArea;

    % --- Accumulate total stress ---
    total_stress_tensor = total_stress_tensor + sigma * area / TotalArea;

end

end


function [vx, vy] = unwrapCellPBC(vx, vy, Lx, Ly)
% PBC-aware (log.txt): a cell straddling the periodic wrap has vertices
% stored far apart in raw coordinates -- polygeom/the edge-length sum
% below on those raw values then blows up for that cell (same root cause
% as TisuePlot.m's "wired" rendering bug). Unwrap every vertex after the
% first relative to the cell's own first vertex (minimum image against
% the true box size Lx,Ly), mirroring ComputeCellColorData.m's
% perCellAreaPerimeter. Harmless no-op for a non-periodic mesh.
if numel(vx) > 1
    dx = vx(2:end) - vx(1); dx = dx - Lx .* round(dx ./ Lx);
    vx(2:end) = vx(1) + dx;
    dy = vy(2:end) - vy(1); dy = dy - Ly .* round(dy ./ Ly);
    vy(2:end) = vy(1) + dy;
end
end


function [ geom, iner, cpmo ] = polygeom( x, y )
%POLYGEOM Geometry of a planar polygon
%
%   POLYGEOM( X, Y ) returns area, X centroid,
%   Y centroid and perimeter for the planar polygon
%   specified by vertices in vectors X and Y.
%
%   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
%   area, centroid, perimeter and area moments of
%   inertia for the polygon.
%   GEOM = [ area   X_cen  Y_cen  perimeter ]
%   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
%     u,v are centroidal axes parallel to x,y axes.
%   CPMO = [ I1     ang1   I2     ang2   J ]
%     I1,I2 are centroidal principal moments about axes
%         at angles ang1,ang2.
%     ang1 and ang2 are in radians.
%     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv

% H.J. Sommer III - 16.12.09 - tested under MATLAB v9.0
%
% sample data
% x = [ 2.000  0.500  4.830  6.330 ]';
% y = [ 4.000  6.598  9.098  6.500 ]';
% 3x5 test rectangle with long axis at 30 degrees
% area=15, x_cen=3.415, y_cen=6.549, perimeter=16
% Ixx=659.561, Iyy=201.173, Ixy=344.117
% Iuu=16.249, Ivv=26.247, Iuv=8.660
% I1=11.249, ang1=30deg, I2=31.247, ang2=120deg, J=42.496
%
% H.J. Sommer III, Ph.D., Professor of Mechanical Engineering, 337 Leonhard Bldg
% The Pennsylvania State University, University Park, PA  16802
% (814)863-8997  FAX (814)865-9693  hjs1-at-psu.edu  www.mne.psu.edu/sommer/

% begin function POLYGEOM

% check if inputs are same size
if ~isequal( size(x), size(y) ),
    error( 'X and Y must be the same size');
end

% temporarily shift data to mean of vertices for improved accuracy
xm = mean(x);
ym = mean(y);
x = x - xm;
y = y - ym;

% summations for CCW boundary
xp = x( [2:end 1] );
yp = y( [2:end 1] );
a = x.*yp - xp.*y;

A = sum( a ) /2;
xc = sum( (x+xp).*a  ) /6/A;
yc = sum( (y+yp).*a  ) /6/A;
Ixx = sum( (y.*y +y.*yp + yp.*yp).*a  ) /12;
Iyy = sum( (x.*x +x.*xp + xp.*xp).*a  ) /12;
Ixy = sum( (x.*yp +2*x.*y +2*xp.*yp + xp.*y).*a  ) /24;

dx = xp - x;
dy = yp - y;
P = sum( sqrt( dx.*dx +dy.*dy ) );

% check for CCW versus CW boundary
if A < 0,
    A = -A;
    Ixx = -Ixx;
    Iyy = -Iyy;
    Ixy = -Ixy;
end

% centroidal moments
Iuu = Ixx - A*yc*yc;
Ivv = Iyy - A*xc*xc;
Iuv = Ixy - A*xc*yc;
J = Iuu + Ivv;

% replace mean of vertices
x_cen = xc + xm;
y_cen = yc + ym;
Ixx = Iuu + A*y_cen*y_cen;
Iyy = Ivv + A*x_cen*x_cen;
Ixy = Iuv + A*x_cen*y_cen;

% principal moments and orientation
I = [ Iuu  -Iuv ;
    -Iuv   Ivv ];
[ eig_vec, eig_val ] = eig(I);
I1 = eig_val(1,1);
I2 = eig_val(2,2);
ang1 = atan2( eig_vec(2,1), eig_vec(1,1) );
ang2 = atan2( eig_vec(2,2), eig_vec(1,2) );

% return values
geom = [ A  x_cen  y_cen  P ];
iner = [ Ixx  Iyy  Ixy  Iuu  Ivv  Iuv ];
cpmo = [ I1  ang1  I2  ang2  J ];

% bottom of polygeom

end


