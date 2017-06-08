function plot_2graphs(data1, data2, block_size, label1, label2)

figure

hold on 
semilogy(abs(2.^data1(:,2)))
semilogy(abs(2.^data2(:,2)))
hold off
title('Change in Pe')

figure;
% Code word length calculation:

hold on 
plot(abs(log2(2.^data1(:,2) ) / block_size )) 
plot(abs(log2(2.^data2(:,2) ) / block_size ))
hold off

legend(label1, label2)
ylim([0 3])
grid on
title('code word length per symbol')

hold off
end