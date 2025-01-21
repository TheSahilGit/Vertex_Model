function [total_stress_tensor, ShearStress_Individual] = calculate_total_stress(Lx,Ly,v,inn,num)


lambda =  1;
A0 = 1;
beta = 0.05;
gamma = 0.06;

para2 = load(strcat("../para2_in.dat"));
Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);

[inside1, inside2, Boundary] = Mesh_Info(Lx,Ly);

TotalArea = 0;
sigma = zeros(2);
total_stress_tensor = zeros(2);

%for ic = 1:Lx*Ly
for ii = 1:length(inside2)
    ic = inside2(ii);

    vx = v(inn(ic,1:num(ic)), 1);
    vy = v(inn(ic,1:num(ic)), 2);
    [ geom, ~, ~ ] = polygeom( vx, vy );
    area = abs(geom(1));

    TotalArea = TotalArea + area;
end

for ii = 1:length(inside2) %1:Lx*Ly

    ic = inside2(ii);

    vx = v(inn(ic,1:num(ic)), 1);
    vy = v(inn(ic,1:num(ic)), 2);

    [ geom, ~, ~ ] = polygeom( vx, vy );
    area = abs(geom(1));
    perimeter = abs(geom(4));

    term1 =  2.0d0 * lambda * (area - A0);
    term2 = (2.0d0 * beta * perimeter + gamma)/(2.0d0*area);

    for j = 1:num(ic)
        jp = j+1;
        if j==num(ic)
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

    sigma(1,1) = term1  + term2 * sigma(1,1);
    sigma(1,2) = term2 * sigma(1,2);
    sigma(2,1) = term2 * sigma(2,1);
    sigma(2,2) = term1 + term2 * sigma(2,2);

    ShearStress_Individual(ii) = sigma(1,2)*area/TotalArea;


end

%for ic = 1:Lx*Ly
for ii = 1:length(inside2)
    ic = inside2(ii);

    vx = v(inn(ic,1:num(ic)), 1);
    vy = v(inn(ic,1:num(ic)), 2);

    [ geom, ~, ~ ] = polygeom( vx, vy );
    area = abs(geom(1));
    total_stress_tensor = total_stress_tensor + sigma * area/TotalArea;

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


