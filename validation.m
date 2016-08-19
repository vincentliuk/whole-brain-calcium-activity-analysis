function [lambdaAmin,lambdaCmin] = validation(lambdaAList,lambdaCList)

% validation part use the first 1000 time series points to 
% validate the hyperparameters lambdaA, lambdaC, by minimizing
% the error ||Y-Y^(lambdaA, lambdaC)||

% call the celegansAnalyze.m to set up

for lambdaA = lambdaAList
    for lambdaC = lambdaCList
        celegansAnalyzeValidation(lambdaA, lambdaC);
    
    end
end

end