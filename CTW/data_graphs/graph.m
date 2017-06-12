
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


fd29 = load('data/Fixed_Depth/Fixed_Depth_NW_004929429.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw29 =  load('data/CTW/CTW_NW_004929429.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw29_comp = load('data/CTW/CTW_NW_004929429.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');


fd30 = load('data/Fixed_Depth/Fixed_Depth_NW_004929430.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw30 =  load('data/CTW/CTW_NW_004929430.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw30_comp = load('data/CTW/CTW_NW_004929430.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');

fd31 = load('data/Fixed_Depth/Fixed_Depth_NW_004929431.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw31 =  load('data/CTW/CTW_NW_004929431.1_S400000_BS200_d16_InvN_CompN.txt','-ASCII');
ctw31_comp = load('data/CTW/CTW_NW_004929431.1_S400000_BS200_d16_InvN_CompY.txt','-ASCII');
%% 28
plot_2graphs(ctw28, fd28, 200, 'CTW', 'Finite Context')
plot_2graphs(ctw28, ctw28_comp, 200, 'CTW', 'CTW - Competitive')
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

python_data = load('data/Sequentially/Sequential_GU170821.1_depth7_InvN.txt', '-ASCII');



Ref_python_data = 2.^load('data/Sequentially/Sequential_GU170821.1_depth7_InvN.txt', '-ASCII');

ref_matlab_data = symb_prob;


diff = ref_matlab_data - Ref_python_data(:,2);

plot(diff)
ylabel('Matlab Prob. - Python Prob. per symbol')
xlabel('Symbol number')

%% Codon vs CTW on coding DNA

codon_data = load('data/Codon_CTW_ENST00000589042.5_MaxS109200_d12_bs300_InvN.txt','-ASCII');
ctw_data = load('data/CTW/CTW_ENST00000589042.5_S109200_BS300_d12_InvN_CompN.txt','-ASCII');

plot_2graphs(codon_data,ctw_data,300,'codon','ctw')

%% Codon vs CTW on Mixed DNA

CTW_codon_DNA = load('data/CTW/CTW_NG_027688.1_S38100_BS300_d12_InvN_CompN.txt','-ASCII');
Codon_codon_DNA = load('data/Codon_CTW_NG_027688.1_MaxS38100_d12_bs300_InvN.txt','-ASCII');

plot_2graphs(Codon_codon_DNA,CTW_codon_DNA, 300, 'Codon','CTW')





