clear; clc; close all;

Lx = 16; % Number of points along x-axis
Ly = 16; % Number of points along y-

israndom = true;
%israndom = false;
%isseed=true;
isseed = false;

%seed = 9324; % remember this seed value
%seed = 9323; % Change this if you want a different random structure.
seed = 9325;


[v,c,p] = MeshGenerator(Lx,Ly,israndom,isseed,seed) ;

v(1,:)=99999;
v(end+64,:) = 0;

[inn,num] = StoreData(Lx,Ly,v,c);
num  = num';


inn(end+64, end+16) = 0; % Extending the dimension of n because in fortran array dimension is static.

Ln = length(num); Li = size(inn,1);
num(Ln+1 : Li) = 0;

para = [Lx;Ly;length(num); length(v(1,:)); length(v(:,1));length(inn(1,:));...
    length(inn(:,1))];

vInitial = v; innInitial = inn; numInitial =  num;

%%

writematrix(v,'v_in.dat','Delimiter',' ');
writematrix(num, 'num_in.dat');
writematrix(inn,'inn_in.dat','Delimiter',' ');
writematrix(para, 'para2_in.dat');



%%



function [borderver] = BorderVertices(Lx,Ly,v,inn,num)
for ii=1:Lx*Ly
    for jj =1:num(ii)
        pk = inn(ii,jj);
        ccn = 0;
        for kk = 1:Lx*Ly
            hh = find(pk==inn(kk, 1:num(kk)),1);
            if hh~=0
                ccn = ccn + 1;
            else
                ccn = ccn + 0;
            end
        end
        cono(ii,jj) = ccn;
    end
end

kk = 1;
for ii = 1:Lx*Ly
    for jj = 1:num(ii)
        if cono(ii,jj) < 3
            cellI(kk) = ii;
            verI(kk)  = jj;
            kk = kk + 1;
        end
    end
end
borderver = [cellI' verI'];

end

function [etaR] = Give_Motility(etasmax, etasmin,Lx,Ly, v,inn,num, Lc)


% intrval = (etasmax - etasmin)/Ly;
%
% for iij = 1 : Ly
%     for i=1:Ly:Lx*Ly-Ly+1
%         etas(i+iij-1, (1:num(i+iij-1))) = etasmax - (iij-1)*intrval;
%     end
% end



for iij = 1 : Ly
    for i = 1:Ly:Lx*Ly-Ly+1
        etas(i+iij-1, (1:num(i+iij-1))) = etasmax*exp(-(iij-1)/Lc);
    end
end

etas(:, end+size(inn,2)-size(etas,2)) = 0;    % Just being lazy. No reason to do this shi!
etas = repmat(etas, 9, 1);
%etas = etas*0;

%%-- Structuring the motility to look like fxx and fyy. So that while using
%%-- it, it does not over add at each vertex.

etaR = zeros(length(v(:,1)),1);

for i = 1:Lx*Ly
    etaR(inn(i, 1:num(i))) = etas(i,1:num(i));
end


end

function [mot] = Give_Motility_new(etasmax, etasmin ,Lx, Ly, v, inn, num, Lc)

for i = 1:Lx*Ly
    for j=1:num(i)
        verry(inn(i,j)) = v(inn(i,j),2);
        verryCount(inn(i,j)) = 1;
    end
end


veryyN = verry(verryCount~=0);
numver = length(veryyN);

mot = zeros(length(v(:,2)),1);

patchwidth = 0.5;

lowpatch = min(veryyN) - patchwidth/2;
highpatch = lowpatch + patchwidth/2;
meanpatch = (highpatch + lowpatch)/2;
for in=1:numver
    for i =1:Lx*Ly
        for j = 1:num(i)
            vy = v(inn(i,j),2);
            if lowpatch <= vy && vy < highpatch          
                mot(inn(i,j)) = etasmax * exp(-vy/Lc);
            end
        end

    end

    lowpatch = highpatch;
    highpatch = lowpatch + patchwidth;
    meanpatch = (highpatch + lowpatch)/2;
end




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
inside2 = setdiff(mainarea,Boundary3);


end



function [v,c,p]=MeshGenerator(Lx,Ly,israndom,isseed,seed)
% Generate the grid of points with hexagonal lattice offset
my = 1:Ly;
mx = 1:Lx;
repeating_sequence = repmat(mx, Ly, 1);
mm = repmat([0 0.5],[Lx,Ly/2])';

x = reshape(repeating_sequence - mm, 1, [])';
y = repmat(my,1,Lx)';


if isseed==true
    rng(seed);
end


ran=rand(1,length(x));
if israndom==true
    x = x(:) ;
    y = y(:) ;

    [inside1, inside2, Boundary] = Mesh_Info(Lx,Ly);

    %        x(inside2') = x(inside2') + 0.5 * (2*ran(inside2')'-1);
    %        y(inside2') = y(inside2') + 0.1 * (2*ran(inside2')'-1);

    x(inside1') = x(inside1') + 0.5 * (2*ran(inside1')'-1);
    y(inside1') = y(inside1') + 0.1 * (2*ran(inside1')'-1);

    x(Boundary') = x(Boundary') + 0.1 * (2*ran(Boundary')'-1);
    y(Boundary') = y(Boundary') + 0.1 * (2*ran(Boundary')'-1);


    %   x = x + 0.5 * (2*ran'-1);
    %   y = y + 0.1 * (2*ran'-1);




    %     for i=2:2:length(x)
    %         x(i)=x(i)+0.5;
    %     end
    %
    %     for i=4:4:length(x)
    %         x(i)=x(i)+0.25;
    %         y(i)=y(i)+0.25;
    %     end

    % y=y(:)+ran(:);
end


% x
% y


xx=[x;x+Lx;x-Lx;x;x;x+Lx;x-Lx;x+Lx;x-Lx];
yy=[y;y;y;y+Ly;y-Ly;y+Ly;y-Ly;y-Ly;y+Ly];
%size(xx)
%[vx,vy]=voronoi(xx,yy);
p=[xx,yy];
%size(p)
[v,c]=voronoin(p);

end

function [inn,num]=StoreData(Lx,Ly,v,c)
for i=1:Lx*Ly
    cellCoord=c{i};
    for j=1:length(cellCoord)
        inn(i,j)=cellCoord(j);
        num(i)=length(cellCoord);
    end
end

%% Restructuring

for i=1:length(num)
    if num(i)~=0
        a=v(inn(i,1:num(i)),1);
        k=find(a==min(a));
        kk=floor(k(1));
        inn(i,1:num(i))=circshift(inn(i,1:num(i)),length(inn(i,1:num(i)))-kk+1);
    end

end

for i=1:length(num)
    if num(i)~=0
        a=v(inn(i,1:num(i)),1);
        if a(2)>a(num(i))
            %disp("yup")
            inn(i,1:num(i))=flip(inn(i,1:num(i)));
            inn(i,1:num(i))=circshift(inn(i,1:num(i)),length(inn(i,1:num(i)))+1);
        end
    end

end

for i=1:length(num)
    if num(i)~=0
        vx=v(inn(i,1:num(i)),1);
        vy=v(inn(i,1:num(i)),2);
        [orient,area]=polyorient(vx,vy);
        if orient
            inn(i,2:num(i))=flip(inn(i,2:num(i)));
            %inn(i,1:num(i))=circshift(inn(i,1:num(i)),length(inn(i,1:num(i)))+1);

        end
    end
end


inn(Lx*Ly+1 :end,:) = 0;
num(Lx*Ly+1 :end) = 0;


% vInitial = v;
% numInitial = num;
% innInitial = inn;
end




%%
function [orient,signed_area] = polyorient(x,y)
%POLYORIENT Orientation of polygon
%   Returns the orientation and signed area of a 2D polygon
%
%   Syntax:
%      [ORIENT,SAREA] = POLYORIENT(X,Y)
%
%   Inputs:
%      X, Y   Vectors with the polygon vertices
%
%   Outputs:
%      ORIENT   Polygon orientation. 1 if the orientation is
%               counter-clockwise (direct), 0 otherwise
%      SAREA    Signed area of the polygon, negative if orientation is
%               not direct
%
%   Examples:
%      x1 = [0 0 1 1]; y1 = [1 2 2 1];
%      x2 = [0 0 1 1]; y2 = [1 0 0 1];
%      x3 = [x1 x2];   y3 = [y1 y2];
%
%      [o1,a1] = polyorient(x1,y1) % 0, -1
%      [o2,a2] = polyorient(x2,y2) % 1,  1
%      [o3,a3] = polyorient(x3,y3) % 0,  0
%
%   MMA 21-4-2006, mma@odyle.net

signed_area=0.5* sum(x.*y([2:end 1]) - y.*x([2:end 1]));
orient = signed_area > 0;
end
