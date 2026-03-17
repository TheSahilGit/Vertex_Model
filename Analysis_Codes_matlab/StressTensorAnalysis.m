clear; clc; close all;


para1 = readtable("../para1_in.dat");
Nt = table2array(para1(7,1));
dt = table2array(para1(8,1));

nrun = 2;

ct = 1;
for it = 1000:1000:7000
    it

    [Lx, Ly, v,inn, num, forces] = LoadData(it, nrun);

    [~, ShearStress_Individual] = calculate_total_stress(Lx,Ly,v,inn,num);

    st(ct) = mean(ShearStress_Individual);
    time(ct) = it*dt;

    ct = ct + 1; 

end

figure()
plot(time, st, '-o')
axis square
xlabel("time")
ylabel("Stress Tensor")
axis square
set(gca, 'FontSize', 28,'LineWidth',2);
set(findall(gca, 'Type', 'Line'), 'LineWidth', 4);

