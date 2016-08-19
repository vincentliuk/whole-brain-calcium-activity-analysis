% set up the lambdaA and lambdaC list 
lambdaAList = [0.00001,0.00005,0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5];
lambdaCList = [0.00001,0.00005,0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5];

% call for the validation function to optimize the best pair of
% lambdaA and lambdaC by minimizing the error between real y and 
% estimated y

[bestLambdaA, bestLambdaC]=validation(lambdaAList, lambdaCList);