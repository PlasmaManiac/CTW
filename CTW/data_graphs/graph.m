
CTW_codon_DNA = load('data/CTW_NG_027688.1_S38100_BS300_d12_InvN_CompN.txt','-ASCII');
Codon_codon_DNA = load('data/Codon_CTW_NG_027688.1_MaxS38100_d12_bs300_InvN.txt','-ASCII');
FD_codon_DNA = load('data/Fixed_Depth_NG_027688.1_S38100_BS300_d12_InvN_CompN.txt','-ASCII');

%%
plot_3graphs(Codon_codon_DNA,CTW_codon_DNA,FD_codon_DNA, 300, 'Codon','CTW','FixedDepth')

%%
plot_2graphs(Codon_codon_DNA,CTW_codon_DNA, 300, 'Codon','CTW')

%% CTW vs FC

fd28 = load('data/Fixed_Depth/Fixed_Depth_NW_004929428.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw28 =  load('data/CTW/CTW_NW_004929428.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw28_comp = load('data/CTW/CTW_NW_004929428.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');

fd28avg = mean(log2(2.^-fd28(:,2)))/200
ctw28avg = mean(log2(2.^-ctw28(:,2)))/200
ctw28_compavg = mean(log2(2.^-ctw28_comp(:,2)))/200 + 0.02


fd29 = load('data/Fixed_Depth/Fixed_Depth_NW_004929429.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw29 =  load('data/CTW/CTW_NW_004929429.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw29_comp = load('data/CTW/CTW_NW_004929429.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');


fd29avg = mean(log2(2.^-fd29(:,2)))/200
ctw29avg = mean(log2(2.^-ctw29(:,2)))/200
ctw29_compavg = mean(log2(2.^-ctw29_comp(:,2)))/200 + 0.02

fd30 = load('data/Fixed_Depth/Fixed_Depth_NW_004929430.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw30 =  load('data/CTW/CTW_NW_004929430.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw30_comp = load('data/CTW/CTW_NW_004929430.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');

fd30avg = mean(log2(2.^-fd30(:,2)))/200 
ctw30avg = mean(log2(2.^-ctw30(:,2)))/200
ctw30_compavg = mean(log2(2.^-ctw30_comp(:,2)))/200 + 0.02


fd31 = load('data/Fixed_Depth/Fixed_Depth_NW_004929431.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw31 =  load('data/CTW/CTW_NW_004929431.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw31_comp = load('data/CTW/CTW_NW_004929431.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');

fd31avg = mean(log2(2.^-fd31(:,2)))/200 
ctw31avg = mean(log2(2.^-ctw31(:,2)))/200
ctw31_compavg = mean(log2(2.^-ctw31_comp(:,2)))/200 + 0.02
%% 28
plot_2graphs(ctw28, fd28, 200, 'CTW', 'Finite Context')
plot_2graphs(ctw28_comp, ctw28,  200, 'CTW - Competitive','CTW')
plot_2graphs(fd28, ctw28_comp, 200, 'Finite Context', 'CTW - Competitive')
%% 29
plot_2graphs(ctw29, fd29, 200, 'CTW', 'Finite Context')
plot_2graphs(ctw29, ctw29_comp, 200, 'CTW', 'CTW - Competitive')
plot_2graphs(fd29, ctw29_comp, 200, 'Finite Context', 'CTW - Competitive')
%% 30
plot_2graphs(ctw30, fd30, 200, 'CTW', 'Finite Context')
plot_2graphs(ctw30, ctw30_comp, 200, 'CTW', 'CTW - Competitive')
plot_2graphs(fd30, ctw30_comp, 200, 'Finite Context', 'CTW - Competitive')
%% 31
plot_2graphs(ctw31, fd31, 200, 'CTW', 'Finite Context')
plot_2graphs(ctw31, ctw31_comp, 200, 'CTW', 'CTW - Competitive')
plot_2graphs(fd31, ctw31_comp, 200, 'Finite Context', 'CTW - Competitive')
%% Comparing values:

python_data_depth7 = load('data/Sequentially/Sequential_GU170821.1_depth7_InvN.txt', '-ASCII');
python_data_depth4 = load('data/Sequentially/Sequential_GU170821.1_depth4_InvN.txt', '-ASCII');

ref_python_data_depth7 = 2.^load('data/Sequentially/Sequential_GU170821.1_depth7_InvN.txt', '-ASCII');
ref_python_data_depth4= 2.^load('data/Sequentially/Sequential_GU170821.1_depth4_InvN.txt', '-ASCII');

ref_matlab_data_depth7 = symb_prob_7;
ref_matlab_data_depth4 = symb_prob_4;


diff7 = ref_matlab_data_depth7 - ref_python_data_depth7(:,2);
diff4 = ref_matlab_data_depth4 - ref_python_data_depth4(:,2);

figure;
plot(diff7)
ylabel('Matlab Prob. - Python Prob. per symbol')
xlabel('Symbol number')
title('Depth 7')

figure;
plot(diff4)
ylabel('Matlab Prob. - Python Prob. per symbol')
xlabel('Symbol number')
title('Depth 4')
%% Codon vs CTW on coding DNA

codon_data = load('data/Codon_CTW_ENST00000589042.5_MaxS109200_d12_bs300_InvN.txt','-ASCII');
ctw_data = load('data/CTW/CTW_ENST00000589042.5_S109200_BS300_d12_InvN_CompN.txt','-ASCII');

plot_2graphs(codon_data,ctw_data,300,'codon','ctw')

codon_data_avg = mean(log2(2.^-codon_data(:,2)))/300 
ctw_data_avg = mean(log2(2.^-ctw_data(:,2)))/300 

%% Codon vs CTW on Mixed DNA

CTW_codon_DNA = load('data/CTW/CTW_NG_027688.1_S38100_BS300_d12_InvN_CompN.txt','-ASCII');
Codon_codon_DNA = load('data/Codon_CTW_NG_027688.1_MaxS38100_d12_bs300_InvN.txt','-ASCII');

plot_2graphs(Codon_codon_DNA,CTW_codon_DNA, 300, 'Codon','CTW')

axis = [1:1:127] * 300;

plot(axis,abs(log2(2.^(Codon_codon_DNA(:,2)/300))))

%% Non - coding DNA

codon_nonC = load('data/Codon_CTW_NW_004929429.1_MaxS30000_d12_bs300_InvN.txt','-ASCII');
ctw_nonC = load('data/CTW_NW_004929429.1_S30000_BS300_d12_InvN_CompN.txt','-ASCII');

codon_nonC_avg = mean(log2(2.^-codon_nonC(:,2)))/200
ctw_non_avg = mean(log2(2.^-ctw_nonC(:,2)))/200

plot_2graphs(codon_nonC,ctw_nonC,300,'Codon', 'CTW')
%% Inverted Complement


inv_28 = load('data/CTW/CTW_NW_004929428.1_S200000_BS200_d16_InvY_CompN.txt','-ASCII');
ctw_28 = load('data/CTW/CTW_NW_004929428.1_S200000_BS200_d16_InvN_CompN.txt','-ASCII');
comp_28 = load('data/CTW/CTW_NW_004929428.1_S200000_BS200_d16_InvN_CompY.txt','-ASCII');
fd_28_inv = load('data/Fixed_Depth/Fixed_Depth_NW_004929428.1_S200000_BS200_d16_InvY_CompY.txt','-ASCII');
comp_inv_28 =  load('data/CTW/CTW_NW_004929428.1_S200000_BS200_d16_InvY_CompY.txt','-ASCII');


avg_inv_28 = mean(log2(2.^-inv_28(:,2)))/400
avg_ctw_28 = mean(log2(2.^-ctw_28(:,2)))/200
comp_inv_28 = mean(log2(2.^-comp_inv_28(:,2)))/400

%% Alpha

ctwA_28 = load('data/CTW/CTW_NW_004929428.1_S400000_BS200_d16_InvN_CompY.txt', '-ASCII');
avg_ctwA_28 =  mean(log2(2.^-ctwA_28(:,2)))/200
plot_2graphs(ctwA_28,fd28,200,'CTW','Finite Context Trees')
