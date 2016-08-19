%function  error = celegansAnalyzeValidation(lambdaA, lambdaC)

%clear all
close all

option = 1;
%% option = 1: apply PLDS algorithm on celegans data and save results
%% option = 2: estimate A & C matrices with PLDS and PCA on first half data
%%
error = 0.0;

% load data
load ../data/Alipasha/3/Supplementary3Data1DeltaF.mat

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
    %y = deltaFToF0; % for active
    %y=Supplementary3Data1DeltaF;
    y=demeanResult_deltaF;
    [p,T]=size(y);
    %m=9;
    m=15;

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
        % initialize C and A matices
        %a = eye(m)*0.5;
        c = U(:,1:m) * sqrt(S(1:m,1:m));
        
        % by using the c and y, reconstruct the x, and then reconsturct ...
        %  A matrix,  calculate x: y= c*x, c'y = c'*c*x, x = c'*c \ c'y
        cty = c'*y;
        ctc = c'*c;
        x = ctc \ cty;
        
        % x(t+1) = A*x(t) 
        xNext = x(:,2:T);
        xCurr = x(:,1:T-1);
        temp1 = xNext*xCurr';
        temp2 = xCurr*xCurr';
        a = temp2 \ temp1 ;

        %a = eye(m)*0.5;
        q=eye(m);
        r=spdiags(ones(p,1),0,p,p);
        Pi=zeros(m,1);
        v=eye(m)*10e-3;
        tol = 10e-3;
        miter = 20;
        
        %[aa,cc,qq,rr,pipi,vv,Sx]=kfs_learn(y,a,c,q,r,Pi,v,tol,miter);
        
        % make estimations & predictions
        penaltyC = linspace(-11,8,19);
        penaltyC = [0 exp(log(4)*penaltyC)];
        %penaltyC = penaltyC([1:2 5:13 15]);
        penaltyA = linspace(-30,-17,14);
        penaltyA = [0 exp(log(4)*penaltyA)];
        penaltyA = penaltyA([1:2 5:13 15]);
        [~,npenal] = size(penaltyC);
        
        accuracy = zeros(1,npenal);
        y_estimate = zeros(p,T);
        x_estimate = zeros(m,T);
        
        minError = 1000000.0;
        penaltySelection = 0.0;
        
        predictStart = 949;
        
        for i = 1:npenal
            lambdaA = penaltyC(i);
            lambdaC = penaltyC(i);
            %lambdaA=100;
            %lambdaC=100;
            % pipip, vvp is used to initialize the x0
            %[aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn(y,a,c,q,r,Pi,v,tol,miter);
            [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y(:,1:predictStart),a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
            
            % initialize x0
            %x_next = mvnrnd(pipip,vvp); % x0
            %x_next = x_next';
            x_next = Sxp(:,predictStart);
            % rrp is new returned R matrix from kfs_learn_p, and rrp is
            % diagonal, sqrtRdiag is taking its sqrt of the diagonal,
            % forming a vector: sqrtRdiag
            sqrtRdiag = sqrt(diag(rrp)); 
%             
%             for j=1:T
%                 
%             end
                
                
            for j = predictStart+1:T
                Q = eye(m);
                ww = mvnrnd(zeros(m,1),Q);
                x_next = aap * x_next + ww';
                %x_next = aap*x_next;
                %x_next = aap * Sxp(:,j-1)+ww'; 
                %x_estimate(:,j) = x_next;
                
                vtmp = normrnd(0,1,p,1);
                vv = sqrtRdiag .* vtmp;
                y_next = ccp * x_next + vv;
                %y_next = ccp * x_next;
                y_estimate(:,j) = y_next;
                
                %y_next = ccp * x_next;
%                 if mod(j,20)==0
%                     figure(j)
%                     plot(y_next,'b');
%                     hold on
%                     plot(y(:,j),'r');
%                     title(['compare real y and estimated y at time point: ',num2str(j)]);
%                     legend('y-recovered','y-origin');
%                     ylabel('signal');
%                     xlabel('neuron#');
%                 end
                    
                
                accuracy(1,i) = accuracy(1,i) + norm((y_next- y(:,j)),2);
                %accuracy(1,i) = accuracy(1,i) + abs(corr(y_next, y(:,j)));
                
            end

            if accuracy(1,i)<=minError
                minError = accuracy(1,i);
                penaltySelection = lambdaC;
            end

            % accuracy(1,i) = sqrt(accuracy(1,i));
            fprintf('lambdaA: %d and lambdaC: %d pair leads to accuracy: %d %n',...
                lambdaA, lambdaC, accuracy(1,i));
        end
        
        fprintf('The selected penalty is: %d', penaltySelection);
        
        % after selecting the penalty, rerun the PLDS
        lambdaA=penaltySelection;
        lambdaC=penaltySelection;
        
        [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y(:,1:predictStart),a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
        
        x_next = Sxp(:,predictStart);
         
        sqrtRdiag = sqrt(diag(rrp)); 
                
        errorPLDS = 0.0;
        
            for j = predictStart+1:T
                Q = eye(m);
                ww = mvnrnd(zeros(m,1),Q);
                x_next = aap * x_next + ww';
                %x_next = aap*x_next;
                %x_next = aap * Sxp(:,j-1)+ww'; 
                %x_estimate(:,j) = x_next;
                
                vtmp = normrnd(0,1,p,1);
                vv = sqrtRdiag .* vtmp;
                y_next = ccp * x_next + vv;
                %y_next = ccp * x_next;
                y_estimate(:,j) = y_next;
                
                %y_next = ccp * x_next;
%                 if mod(j,20)==0
%                     figure(j)
%                     plot(y_next,'b');
%                     hold on
%                     plot(y(:,j),'r');
%                     title(['compare real y and estimated y at time point: ',num2str(j)]);
%                     legend('y-recovered','y-origin');
%                     ylabel('signal');
%                     xlabel('neuron#');
%                 end
                    
                
                errorPLDS = errorPLDS + norm((y_next- y(:,j)),2);
                %accuracy(1,i) = accuracy(1,i) + abs(corr(y_next, y(:,j)));
                
            end
        
        
        % draw sample prediction traces for selected lambda (PLDS)
        for neuron = 1:15:99
            figure(neuron)
            plot(predictStart+1:T,y(neuron,predictStart+1:T),'r');
            hold on
            plot(predictStart+1:T,y_estimate(neuron,predictStart+1:T),'b');
            title(['Prediction:compare real y and estimated y traces for neuron#: ',num2str(neuron)]);
            legend('y-origin','y-recovered');
            ylabel('signal');
            xlabel('time (0.2s)');
        end
        
        % apply pca parameters to the prediction.
        % the pca parameters are 'a' (A matrix), and 'c' (C matrix)
        accuracyPCA = 0.0;
        x_PCAnext = x(:,predictStart);
        y_PCAestimate = zeros(p,T);
        for j = predictStart+1:T
                Q = eye(m);
                wwPCA = mvnrnd(zeros(m,1),Q);
                x_PCAnext = aap * x_PCAnext+wwPCA';
                
                vtmpPCA = normrnd(0,1,p,1);
                y_PCAnext = ccp * x_PCAnext+vtmpPCA;
                y_PCAestimate(:,j) = y_PCAnext;
            
                accuracyPCA = accuracyPCA + norm((y_PCAnext- y(:,j)),2);
                %accuracy(1,i) = accuracy(1,i) + abs(corr(y_next, y(:,j)));
                
        end
        % draw samples prediction from PCA estimate
        for neuron = 1:15:99
            figure(neuron+100)
            plot(predictStart+1:T,y(neuron,predictStart+1:T),'r');
            hold on
            plot(predictStart+1:T,y_PCAestimate(neuron,predictStart+1:T),'b');
            title(['PCA Prediction:compare real y and estimated y traces for neuron#: ',num2str(neuron)]);
            legend('y-origin','y-recovered');
            ylabel('signal');
            xlabel('time (0.2s)');
        end
   
        
        
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

% calculate the error between the real y and estimated y by 
% using the estimated C and A matrices 

%end

