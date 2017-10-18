function plot_2graphs(data1, data2, block_size, label1, label2)

%figure

hold on 
%semilogy(abs(2.^data1(:,2)))
%semilogy(abs(2.^data2(:,2)))
hold off
title('Change in Pe')

axis = [1:1:length(data1(:,2))] * block_size;

figure;
% Code word length calculation:

plot(axis,abs(log2(2.^data1(:,2) ) / (block_size  ) ), 'linewidth',1.5) 
hold on 
plot(axis,abs(log2(2.^data2(:,2) ) / (block_size) ), 'linewidth',0.5)
hold off

legend(label1, label2)
grid on
title('Average Code Word Length per Symbol')
xlabel('Symbol')
ylabel('Average Code Word Length [bits]')


figure;
plot(axis,data2(:,1))
title('Best Depth per block')
xlabel('Depth')
ylabel('Average Code Word Length [bits]')




hold off
end