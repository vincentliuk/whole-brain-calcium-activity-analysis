%clear all
close all

load ../data/rescaledata.mat

y = extSignals_rescale;
% we found that extSignals_norm_ValidOnly data does not have downward
% spikes
y = extSignals_norm_ValidOnly;

for i = 1 : 100
    figure(1)
    hold on
    plot(y(i,:));
    xlabel('n (frame)')
    ylabel('signal')
end

for i = 1 : 10
    figure(2)
    hold on
    plot(y(i,:));
    xlabel('n (frame)')
    ylabel('signal')
end