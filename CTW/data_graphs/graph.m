
rand16 = load('ctw_random_16.txt','-ASCII');

%%
fd_inv = load('data/Fixed_Depth_L_S99900_BS300_d12_InvY_CompN.txt','-ASCII');

block_size = 300;

%%
codon12_inv = load('data/Codon_CTW_L_MaxS99900_d12_bs300_InvY.txt','-ASCII');

%%
plot_graphs(ctw, 300)
plot_graphs(fd, block_size)
plot_2graphs(ctw, fd , block_size)

%%
rand_noInv = load('ctw_random_blockNoInverse_noA.txt','-ASCII');

%%
plot_3graphs(codon12_inv, ctw_incv, fd_inv, 600, 'codon', 'CTW', 'Fixed Depth')