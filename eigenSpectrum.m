%close all
%clear all

%load ../data/SVDresult.mat

figure(1)

eigen = diag(S);
eigen = eigen.^2;

plot(eigen)
xlim([0 50])
xlabel('n')
ylabel('\lambda_n')
set(gca,'yscale','log')