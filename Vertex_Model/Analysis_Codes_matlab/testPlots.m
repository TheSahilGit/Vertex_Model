clear; clc ;
%close all;


para2 = load("para2_in.dat");
etas = load("motility_in.dat");



%%

Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);



for it = 99999000  
    it
    im = floor(it/1000)  % 100 is the data writing interval. Change it as per the Fortran code.
    fname_inn = sprintf('data/inn_%08d.dat', it);
    fid = fopen(fname_inn);
    dum1 = fread(fid,1,'float32');
    inn = fread(fid,inndim1*inndim2,'int');
    inn = reshape(inn, [inndim1,inndim2]);
    inn = inn';

    fname_num = sprintf('data/num_%08d.dat', it);
    fid = fopen(fname_num);
    dum2 = fread(fid,1,'float32');
    num = fread(fid,numdim,'int');
    num(Lx*Ly+1:end) = [];

    fname_v = sprintf('data/v_%08d.dat', it);
    fid = fopen(fname_v);
    dum3 = fread(fid,1,'float32');
    v = fread(fid,vdim1*vdim2,'float64');
    v = reshape(v,[vdim1,vdim2]);
    v = v';

    TisuePlot(Lx,Ly,v,inn,num,etas)

 
end

