function plot_graphs(data)
figure
subplot(2,1,1)
plot(data(:,1))
ylim([0 16]);
title('Best Depth')

subplot(2,1,2)
semilogy(abs(2.^data(:,2)))
title('Change in Pe')

figure;
% Code word length calculation:
plot(abs(log2(2.^data(:,2) ) / 400 )) 
ylim([0 3])
grid on
title('code word length per symbol')
end