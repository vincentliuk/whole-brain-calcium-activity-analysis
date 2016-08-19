close all
% calculate x: y= c*x, c'y = c'*c*x, x = c'*c \ c'y
cty = ccp'*y;
ctc = ccp'*ccp;
x2 = ctc \ cty;
x2 = x2';
Ynew = ccp * x2';

imagesc(Ynew);
caxis([-1 5]); % set color bar range from -1~5
title('Ynew');
xlabel('time frames');
ylabel('neuron No.');

diff = Ynew - y;
diff_relative = abs(diff) ./ y;

figure(1)
imagesc(diff);

diff_norm = norm(diff,2);

y_norm = norm(y,2);


for i =2 :6
    figure(i)
    plot(x2(:,i));
end

% plot some example trajectory of y
for i = 1:20;
    figure(10+i)
    plot(y(i,:),'r');
    hold on;
    plot(Ynew(i,:),'b');
    legend('y-origin','y-recovered');
    filename = strcat('y_compare_plot',num2str(i),'.jpg');
    ylabel('y');
    xlabel('time');
    print(filename,'-djpeg')
end
