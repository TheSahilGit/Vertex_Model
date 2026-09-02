clear; clc; close all;


para1 = readtable("../para_Simulation.dat");
Nt = table2array(para1(7,1));
dt = table2array(para1(8,1));

nrun = 2;

ct = 1;
for it = 100:100:10000
    it

    [Lx, Ly, v,inn, num, forces] = LoadData(it, nrun);
    
    radius = 10; 
    [~, ShearStress_Individual] = Calculate_Total_Stress(Lx,Ly,v,inn,num, radius);

    st(ct) = mean(ShearStress_Individual);
    time(ct) = it*dt;

    ct = ct + 1; 

end

%writematrix([time' st'], 'st.dat')

figure()
plot(time, st, '-o')
axis square
xlabel("time")
ylabel("Stress Tensor")
axis square
set(gca, 'FontSize', 28,'LineWidth',2);
set(findall(gca, 'Type', 'Line'), 'LineWidth', 4);

