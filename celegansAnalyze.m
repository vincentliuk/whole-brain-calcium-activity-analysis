%clear all
close all

option = 1;
%% option = 1: apply PLDS algorithm on motor data and save results
%% option = 2: estimate A & C matrices with PLDS and PCA on first half data
%%

load ../data/Alipasha/3/Supplementary3Data1DeltaF.mat % load data

if option == 1
    %   analyze real data with KFS and penalized KFS model
    %   KIRBY 21 Motor Network Analysis

    %   import data
  %  yy=[];
    %   datapath='./data/motorData/';
    %   use fsl smoothed data
    
  %  datapath = './data/';
  %  yy(1,:,:)=csvread([datapath,'KKI2009-01.csv'],1,0);
  %  yy(2,:,:)=csvread([datapath,'KKI2009-25.csv'],1,0);
  %  yy(3,:,:)=csvread([datapath,'KKI2009-04.csv'],1,0);
  %  yy(4,:,:)=csvread([datapath,'KKI2009-11.csv'],1,0);
    %   yy(5,:,:)=csvread([datapath,'KKI2009-02.csv'],1,0);
    %   yy(6,:,:)=csvread([datapath,'KKI2009-37.csv'],1,0);
    %   relavant dimensions
 %   y = extSignals_rescale;
    % y = demeanRest; % for resting 
    y = demeanResult_deltaFToF0; % for active
    [p,T]=size(y);
    %m=9;
    m=20;

    %   results containers
    kfsA={};
    kfsC={};
    pkfsA={};
    pkfsC={};

    %   analyze
%    for i = 1:4
%        y = squeeze(yy(i,:,:))';
        i = 1;
        [U,S,V] = svd(y,'econ');
        a = eye(m)*0.5;
        c = U(:,1:m) * sqrt(S(1:m,1:m));
        q=eye(m);
        r=spdiags(ones(p,1),0,p,p);
        Pi=zeros(m,1);
        v=eye(m)*10e-3;
        tol = 10e-3;
        miter = 20;


        %lambdaA = 10000;
        %lambdaC = 0.1;
        %lambdaA = 20000;
        %lambdaC = 0.25;
        lambdaA = 0;
        lambdaC = 0;
        %   lambdaA: below 2.5e-5
        %lambdaA = 0.0000258;
        lambdaA = 0.000000001;
        %   lambdaC: from 0/1e-9 to 1e-6
        %lambdaC = 0.000000001;
        lambdaC = 0.000000002;

        %[aa,cc,qq,rr,pipi,vv,Sx]=kfs_learn(y,a,c,q,r,Pi,v,tol,miter);
        [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
        %kfsA{i}=aa;
        %kfsC{i}=cc;
        pkfsA{i}=aap;
        pkfsC{i}=ccp;
%    end
    %   save('./results/motor/motor_kfs.mat','kfsA','kfsC')
    %   save('./results/motor/motor_pkfs.mat','pkfsA','pkfsC')
    %save('./results/motor/smoothed/motor_kfs.mat','kfsA','kfsC')

    %save('./results/motor/smoothed/motor_pkfs.mat','pkfsA','pkfsC')
    %save('./results/demean_active_pkfs_5.mat');
    %save('./results/celegans_pkfs_5.mat','pkfsA','pkfsC')
end
if option == 2

    %   import data
  %  yy=[];
    %   use fsl smoothed data
  %  datapath = './data/motorData/smoothed/';
  %  yy(1,:,:)=csvread([datapath,'KKI2009-01.csv'],1,0);
    %   relavant dimensions
    yy = extSignals_rescale;
    [p,T]=size(y);
    m=5;
    % m is the number of latent variable dimension, i.e., dimension of x
    train_t = T/2;
    
    Y = yy;
    y = Y(:,1:train_t);
    [U,S,V] = svd(y,'econ');
    a = eye(m);
    c = U(:,1:m) * S(1:m,1:m);
    q=eye(m);
    r=spdiags(ones(p,1),0,p,p);
    Pi=zeros(m,1);
    v=eye(m)*10e-3;
    tol = 10e-3;
    miter = 20;

    lambdaA = 0.1;
    lambdaC = 0.1;   
    
    [aa,cc,qq,rr,pipi,vv,Sx]=kfs_learn(y,a,c,q,r,Pi,v,tol,miter);
    [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
    
    % pca algorithm
    % this part will be done in R
    
    %save('./results/celegans_kfs_pkfs_26_pred.mat','Y','y','train_t','U','V','S','aa','cc','aap','ccp','Sx','Sxp')
    
end

