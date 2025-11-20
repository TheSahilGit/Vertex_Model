
function [Lx, Ly, v,inn,num, forces, if_alive, T2_time] = LoadData(it, nrun)


para2 = load("../para2_in.dat");
etas = load("../motility_in.dat");




%%

Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);



if nrun==1
    fname_inn = sprintf('../data/inn_%08d.dat', it);
    fname_num = sprintf('../data/num_%08d.dat', it);
    fname_v = sprintf('../data/v_%08d.dat', it);
    fname_force = sprintf('../data/force_%08d.dat', it);
elseif nrun==2
    fname_inn = sprintf('../data/nrun2_inn_%08d.dat', it);
    fname_num = sprintf('../data/nrun2_num_%08d.dat', it);
    fname_v = sprintf('../data/nrun2_v_%08d.dat', it);
    fname_force = sprintf('../data/nrun2_force_%08d.dat', it);
end

fid = fopen(fname_inn);
dum1 = fread(fid,1,'float32');
inn = fread(fid,inndim1*inndim2,'int');
inn = reshape(inn, [inndim1,inndim2]);
inn = inn';

fid = fopen(fname_num);
dum2 = fread(fid,1,'float32');
num_1 = fread(fid,numdim*3,'int');
num_1 = reshape(num_1,[3, numdim]);
num = num_1(1,:);
if_alive = num_1(2,:);
T2_time = num_1(3,:);
%num(Lx*Ly+1:end) = [];

%

fid = fopen(fname_v);
dum3 = fread(fid,1,'float32');
v = fread(fid,vdim1*vdim2,'float64');
v = reshape(v,[vdim1,vdim2]);
v = v';



% fid = fopen(fname_force);
% dum4 = fread(fid,1,'float32');
% forces = fread(fid,vdim1*vdim2,'float64');
% forces = reshape(forces,[vdim1,vdim2]);
% forces = forces';


fid = fopen(fname_force);
dum4 = fread(fid,1,'float32');
forces = fread(fid,8*vdim2,'float64');
forces = reshape(forces,[8,vdim2]);
forces = forces';






end
