load("msd.dat");

figure()

loglog(msd(:,1), msd(:,2), '-o')
hold on
loglog(msd(:,1), msd(:,3), '-o')
hold on
loglog(msd(:,1), msd(:,4).*1e2, '-o')

