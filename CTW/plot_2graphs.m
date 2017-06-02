function plot_2graphs(data1, data2)

figure

subplot(2,1,1)
hold on 
plot(data1(:,1))
plot(data2(:,1))
hold off
ylim([0 16]);
title('Best Depth')

subplot(2,1,2)

hold on 
semilogy(abs(2.^data1(:,2)))
semilogy(abs(2.^data2(:,2)))
hold off
title('Change in Pe')

figure;
% Code word length calculation:

hold on 
plot(abs(log2(2.^data1(:,2) ) / 400 )) 
plot(abs(log2(2.^data2(:,2) ) / 400 ))
hold off

ylim([0 3])
grid on
title('code word length per symbol')

hold off
end