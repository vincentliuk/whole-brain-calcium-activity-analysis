clear all
close all


load ../data/active.mat
load ../data/extSignals_Recon3D_20140311_fish2_20Hz_stim1min_20X_05NA.mat
% here is the preprocssing of the neural data
% substract mean:
% for each neuron, we substracted the mean value averaged over the
% measurement time period. 


active = extSignals_norm_ValidOnly(:,301:end);
rest = extSignals_norm_ValidOnly(:,1:300);
restMean = mean(rest,2);

for i = 1 : 10
    figure(1)
    hold on
    plot(rest(i,:))
    xlabel('n (frame)')
    ylabel('signal')
end

ntime = size(rest,2);
demeanRest = rest - repmat(restMean,1,ntime);


meanval = mean(active,2);

for i = 1 : 10
    figure(2)
    hold on
    plot(active(i,:));
    xlabel('n (frame)')
    ylabel('signal')
end

numNeuron = size(active,1);
ntime = size(active,2);
demean = active - repmat(meanval,1,ntime);

for i = 1 : 10
    figure(3)
    hold on
    plot(demean(i,:));
    xlabel('n (frame)')
    ylabel('signal')
end

% after the mean, we want to:
% at each time point, we substrate the mean average over the all the
% neurons at this time point and then divided by the sqrt root of variance

mean_over_neuron = mean(active);
demean_over_neuron = active - repmat(mean_over_neuron,numNeuron,1);

for i = 1 : 10
    figure(3)
    hold on
    plot(demean_over_neuron(i,:));
    xlabel('n (frame)')
    ylabel('signal')
end

std_over_neuron = std(active);
rescale_over_neuron = demean_over_neuron./repmat(std_over_neuron,numNeuron,1);

for i = 1 : 10
    figure(4)
    hold on
    plot(rescale_over_neuron(i,:));
    xlabel('n (frame)')
    ylabel('signal')
end
