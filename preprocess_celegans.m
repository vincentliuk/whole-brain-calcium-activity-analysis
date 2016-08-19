%clear all
%close all

load ../data/Alipasha/3/Supplementary3Data1DeltaF
load ../data/Alipasha/3/Supplementary3Data1FLraw

% here is the preprocssing of the neural data
% substract mean:
% for each neuron, we substracted the mean value averaged over the
% measurement time period. 


%active = extSignals_norm_ValidOnly(:,301:end);
%rest = extSignals_norm_ValidOnly(:,1:300);
%restMean = mean(rest,2);

% for i = 1 : 10
%     figure(1)
%     hold on 
%     plot(rest(i,:))
%     xlabel('n (frame)')
%     ylabel('signal')
% end

meanData_deltaF = mean(Supplementary3Data1DeltaF,2);
ntime_deltaF = size(Supplementary3Data1DeltaF,2);
demeanResult_deltaF = Supplementary3Data1DeltaF - repmat(meanData_deltaF,1,ntime_deltaF);
std_over_time_deltaF = std((Supplementary3Data1DeltaF'));

% divide each deltaF by its standard deviation 
for i=1:size(Supplementary3Data1DeltaF,1)
    demeanResult_deltaF(i,:) = demeanResult_deltaF(i,:)./std_over_time_deltaF(i);
end


figure(1)
hold on
imagesc(Supplementary3Data1DeltaF);
title('3\_SupplementaryData1\_S1\_DeltaF');
figure(2)
hold on
imagesc(demeanResult_deltaF);
title('3\_SupplementaryData1\_S1\_DeltaF\_demeanResult');

figure(3)
hold on
plot(meanData_deltaF);
title('mean of each neuron traces')
xlabel('neuron #')
ylabel('deltaF mean')


figure(4)
hold on
plot(std_over_time_deltaF);
title('standard deviation among all neurons over time series')
xlabel('neuron (#)')
ylabel('deltaF mean')

figure(5)
hold on
imagesc(Supplementary3Data1FLraw);
title('3\_SupplementaryData1\_FLraw');


% calculate deltaF/F0 = 100*DeltaF/(FLraw-DeltaF)
deltaFToF0 = 100.*(Supplementary3Data1DeltaF)./(Supplementary3Data1FLraw-Supplementary3Data1DeltaF);
figure(6)
hold on
imagesc(deltaFToF0);
title('3\_SupplementaryData1\_deltaF/F0');


meanData_deltaFToF0 = mean(deltaFToF0,2);
ntime_deltaFToF0 = size(deltaFToF0,2);
demeanResult_deltaFToF0 = deltaFToF0 - repmat(meanData_deltaFToF0,1,ntime_deltaFToF0);

% meanval = mean(active,2);
% 
% for i = 1 : 10
%     figure(2)
%     hold on
%     plot(active(i,:));
%     xlabel('n (frame)')
%     ylabel('signal')
% end
% 
% numNeuron = size(active,1);
% ntime = size(active,2);
% demean = active - repmat(meanval,1,ntime);
% 
% for i = 1 : 10
%     figure(3)
%     hold on
%     plot(demean(i,:));
%     xlabel('n (frame)')
%     ylabel('signal')
% end