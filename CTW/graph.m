
ctw = load('CTW.txt','-ASCII');
markov = load('markov.txt','-ASCII');

figure
subplot(2,1,1)
plot(ctw(:,1))
ylim([0 18]);
title('Best Depth - CTW')

subplot(2,1,2)
plot(abs(ctw(:,2)))
title('Change in Pw - CTW')

figure;
% Code word length calculation:
plot(abs(log2(2.^ctw(:,2) ) / 200 ))
ylim([0 3])
grid on
title('code word length - CTW')

%% Fixed Depth
figure
subplot(2,1,1)
plot(markov(:,1))
ylim([0 16]);
title('Best Depth - FD')

subplot(2,1,2)
plot(abs(markov(:,2)))
title('Change in Pe - FD')

figure;
% Code word length calculation:
plot(abs(log2(2.^markov(:,2) ) / 200 )) 
ylim([0 3])
grid on
title('code word length- FD')
