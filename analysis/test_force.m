clear; clc; close all;



nrun = 2;

ct = 0;
for it = 10000:10000:1000000
    ct = ct + 1;
    [Lx, Ly, v,inn,num, forces] = LoadData(it, nrun);


    idx = find(inn(:,1) == 0, 1, 'first');


    fxx = forces(:,1);
    fyy = forces(:,2);
    fxx_ran = forces(:,3);
    fyy_ran = forces(:,4);
    fxx_active_contr = forces(:,5);
    fyy_active_contr = forces(:,6);
    fxx_ABP = forces(:,7);
    fyy_ABP = forces(:,8);


    fxtotal =  fxx + fxx_ran + fxx_active_contr + fxx_ABP;
    fytotal =  fyy + fyy_ran + fyy_active_contr + fyy_ABP;


    ftotal(ct) = max(sqrt(fxtotal.^2 + fytotal.^2));
    f_contr_total = sqrt(fxx_active_contr.^2 + fyy_active_contr.^2);
    f_det_total = sqrt(fxx.^2 + fyy.^2);

end


% disp('max')
% max(ftotal)
% disp('max active contr')
% max(f_contr_total)
% disp('max det ')
% max(f_det_total)
% 

figure()
semilogy(ftotal)
axis square



