% author: shaojie chen (schen89@jhu.edu)
% plds simulation 2

% use simdata_gen.m to generate data first
load('p-300-d-10-T-100-data-regA-smoothsparseC.mat');

penaltyC = linspace(-11,2,14);
penaltyC = [0 exp(log(4)*penaltyC)];
penaltyA = linspace(-30,-17,14);
penaltyA = [0 exp(log(4)*penaltyA)];
[~,npenal] = size(penaltyC);

cormat={};

for i = 1:npenal
    lambdaA = penaltyA(i);
    lambdaC = penaltyC(i);
    
    [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
    cormat{i,1}=corr(A,aap);
    cormat{i,2}=corr(C,ccp);
    
    if i == 1
        aa = aap;
        cc = ccp;
    end
    
    if i == 7
        paa = aap;
        pcc = ccp;
    end
end

save(['./p-',num2str(p),'-d-',num2str(d),'-T-',num2str(T),'-sim2-result.mat'],'p','d','T','penaltyA','penaltyC','cormat','A','C','aa','cc','paa','pcc');

% next step: use the simulation1.R for further analysis and plots
