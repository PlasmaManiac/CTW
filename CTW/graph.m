close all;

ctw_random_1A = load('update/ctw_random_block_-1_A.txt','-ASCII');

fd_random_1A = load('update/fd_random_block_-1_A.txt','-ASCII');
%%


plot_graphs(ctw_40000_31_A1)
%plot_graphs(fd_40000_31_A1)

plot_graphs(ctw_40000_31_A)
%plot_graphs(fd_40000_31_A)

plot_graphs(ctw_40000_31_1)
%plot_graphs(fd_40000_31_1)

%% Max Depth = 6

plot_graphs(ctw_MD6)
plot_graphs(markov_MD6)

%% Max Depth = 12
plot_graphs(ctw_MD12)
plot_graphs(markov_MD12)
%% Repeating A
ctw_repeat = load('repeating_A_CTW.txt', '-ASCII');
fd_repeat  = load('repeating_A_FD.txt', '-ASCII');

figure;
semilogy(abs(ctw_repeat(:,2)))
grid on
title('Change in Pw - CTW')

figure;
semilogy(abs(log2(2.^ctw_repeat(:,2) ) / 200 ))
grid on
title('code word length - CTW')


figure;
semilogy(abs(fd_repeat(:,2)))
grid on
title('Change in Pw - FD')


figure;
semilogy(abs(log2(2.^fd_repeat(:,2) ) / 200 ))
grid on
title('code word length - FD')


%% Repeating ACGTT

ctw_repeat = load('repeating_ACGTT_CTW.txt', '-ASCII');
fd_repeat  = load('repeating_ACGTT_FD.txt', '-ASCII');

figure;
semilogy(abs(ctw_repeat(:,2)))
grid on
title('Change in Pw - CTW')

figure;
semilogy(abs(log2(2.^ctw_repeat(:,2) ) / 200 ))
grid on
title('code word length - CTW')

figure;
semilogy(abs(fd_repeat(:,2)))
grid on
title('Change in Pw - FD')

figure;
semilogy(abs(log2(2.^fd_repeat(:,2) ) / 200 ))
grid on
title('code word length - FD')

%% Random Block
ctw_random = load('randomBlock_CTW.txt', '-ASCII');
fd_random  = load('randomBlock_FD.txt', '-ASCII');

figure;
plot(abs(ctw_random(:,2)))
grid on
title('Change in Pw - CTW')
xlim([0 30])

figure;
plot(abs(log2(2.^ctw_random(:,2) ) / 200 ))
grid on
title('code word length - CTW')
xlim([0 30])

figure;
plot(abs(fd_random(:,2)))
grid on
title('Change in Pw - FD')
xlim([0 30])


%% Random Blocks
ctw_randoms = load('randomBlocks_CTW.txt', '-ASCII');
fd_randoms  = load('randomBlocks_FD.txt', '-ASCII');

figure;
plot(abs(ctw_randoms(:,2)))
grid on
title('Change in Pw - CTW')

figure;
plot(abs(log2(2.^ctw_randoms(:,2) ) / 200 ))
grid on
title('code word length - CTW')

figure;
plot(abs(fd_randoms(:,2)))
grid on
title('Change in Pw - FD')

figure;
plot(abs(log2(2.^fd_randoms(:,2) ) / 200 ))
grid on
title('code word length - FD')
%%
figure;
plot(abs(ctw0001(:,2)))
hold on;
plot(abs(ctw001(:,2)))
plot(abs(ctw01(:,2)))
plot(abs(ctw05(:,2)))

grid on
title('Change in Pw - CTW')
legend('.0001','.001','.01','.05')
%%
figure;
plot(abs(log2(2.^ctw0001(:,2) ) / 200 ))
hold on;
plot(abs(log2(2.^ctw001(:,2) ) / 200 ))
plot(abs(log2(2.^ctw01(:,2) ) / 200 ))
plot(abs(log2(2.^ctw05(:,2) ) / 200 ))


ylim([0 1.5])
grid on
title('code word length - CTW')
legend('.0001','.001','.01','.05')
%%
figure
subplot(2,1,1)
plot(ctw(:,1))
ylim([0 18]);
title('Best Depth - CTW')
grid on;
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
