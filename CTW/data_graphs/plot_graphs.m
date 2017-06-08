function plot_graphs(data, block_size)
%  figure
%  subplot(2,1,1)
%  plot(data(:,1))
%  ylim([0 16]);
%  title('Best Depth')
% 
% subplot(2,1,2)
figure;
semilogy(abs(2.^data(:,2)))
title('Change in Pw')
grid on;

figure;
% Code word length calculation:
plot(abs(log2(2.^data(:,2)) / block_size)) 
ylim([0 2.5])
grid on
title('code word length per symbol')
end