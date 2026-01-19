clear; clc; close all;



para1 = readtable("../para1_in.dat");
dt = table2array(para1(8,1));
nrun = 2;


itInit = 1000; 
interval = 1000; 
itfinal = 400000; 

ac = 0.7; 
[time, Qt] = CalculateQt(itInit,interval, itfinal, dt, nrun, ac);

%writematrix([time' Qt'], 'Q-time_lowmot.dat');



%%
figure("Position",[100 100 800 800])
semilogx(time, Qt,'o', "LineWidth",3, 'MarkerSize',20);
%hold on; 
xlabel("Time");
ylabel("Q(t)")
set(gca, "FontSize", 32)
axis square
legend()





function [time, Qt] = CalculateQt(itInit,interval, itfinal, dt, nrun, ac)

%ac = 1.414; 
time = 0;
ct = 1;
for it = itInit:interval:itfinal
    

    [Lx, Ly, v, inn, num, ~, ~] = LoadData(it, nrun);


    for ic = 1:Lx*Ly
        vx = v(inn(ic, 1:num(ic)),1);
        vy = v(inn(ic, 1:num(ic)),2);

        vcmX(ic) = mean(vx);
        vcmY(ic) = mean(vy);

    end

    if it == itInit
        vcmXIn = vcmX;
        vcmYIn = vcmY;
    end

    W = 0; 
    for ic = 1:Lx*Ly
        dis = sqrt((vcmX(ic) - vcmXIn(ic)).^2 + (vcmY(ic) - vcmYIn(ic)).^2);

        if ac - dis  > 0
            W = W + 1;
        end
    end

    Qt(ct) = W/(Lx*Ly);


    time(ct) = it*dt;

    ct = ct + 1

    clear v inn num
end

end