function myplot()

load ./results/iteration.dat;
iteration;

%%%% energy
%plot(iteration(:,3));

figure1 = figure('color', 'white');
axes1 = axes('Parent', figure1, 'FontName','Times New Roman',...
	 'FontSize',30);
set(gca,'linewidth',3)      %| 设置图形外边框的线宽
set(gcf, 'unit', 'normalized', 'position', [0.05, 0.05, 0.85, 0.85])
box on

%%%% energy difference
semilogy(abs(iteration(:,4)), 'LineStyle', '-', 'LineWidth', 2, ...
	'color', [0 0 0]);
hold on 

ylimLeft = 1;
ylimRight = 5.0e-07;
%ylim([ylimLeft, ylimRight])

%%% ITER 0 : Radius = 3.650000000000000e+00, energy = -2.949848340868652e+00
%%% ITER 1 : Radius = 4.811290147411177e+00, energy = -2.986835998703285e+00, dtol_energy = -3.698765783463287e-02, 
%%% ITER 2 : Radius = 4.352960596564764e+00, energy = -2.991832958819238e+00, dtol_energy = -4.996960115953186e-03, 
%%% ITER 3 : Radius = 4.462402959533113e+00, energy = -2.992482580902515e+00, dtol_energy = -6.496220832774569e-04, 
%%% ITER 4 : Radius = 4.427244908949050e+00, energy = -2.992426641390829e+00, dtol_energy = 5.593951168592071e-05, 
%%% ITER 5 : Radius = 4.438065810468289e+00, energy = -2.992371444291748e+00, dtol_energy = 5.519709908163861e-05, 
%%% ITER 6 : Radius = 4.434223940194329e+00, energy = -2.992448608108361e+00, dtol_energy = -7.716381661371230e-05, 
%%% ITER 7 : Radius = 4.436024745554335e+00, energy = -2.992439908185184e+00, dtol_energy = 8.699923177424296e-06, 
%%% ITER 8 : Radius = 4.435214141208949e+00, energy = -2.992452212393655e+00, dtol_energy = -1.230420847075919e-05, 
%%% ITER 9 : Radius = 4.435574689038679e+00, energy = -2.992458565863010e+00, dtol_energy = -6.353469355069308e-06, 
%%% ITER 10 : Radius = 4.435621488760628e+00, energy = -2.992449293589280e+00, dtol_energy = 9.272273730154268e-06, 
%%% ITER 11 : Radius = 4.435475019647649e+00, energy = -2.992445956436546e+00, dtol_energy = 3.337152733706716e-06, 
%%% ITER 12 : Radius = 4.435428206848609e+00, energy = -2.992445839048598e+00, dtol_energy = 1.173879482507800e-07, 


end
