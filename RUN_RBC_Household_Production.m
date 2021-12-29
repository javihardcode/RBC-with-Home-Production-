%% Run: RBC With Household Production: 
clear all; close all; 

addpath '/Applications/Dynare/4.6.3/matlab' ; 
close all; 

dynare MOD_RBC_Household_Production.mod 



%% Here we run experiments and figures in Dynare: 
% variables to represent: Y k hm hf lm lf H 

% Shock to the Market good: 

figure; 
subplot(3,2,1); 
plot(Y_ee, 'k');
hline = refline(0, 0);
title('GDP'); 

subplot(3,2,2); 
plot(k_ee, 'k');
hline = refline(0, 0);
title('Capital'); 

subplot(3,2,3); 
plot(hm_ee, 'k');
hline = refline(0, 0);
title('Male Market Labor Supply'); 

subplot(3,2,4); 
plot(hf_ee, 'k'); 
hline = refline(0, 0);
title('Female Market Labor Supply');


subplot(3,2,5); 
plot(lm_ee, 'k');
hline = refline(0, 0);
title('Male Home Labor Supply'); 

subplot(3,2,6); 
plot(lf_ee, 'k'); 
hline = refline(0, 0);
title('Female home Labor Supply');



close all; 

% Shock to the Market good: 
figure; 
subplot(3,2,1); 
plot(Y_uu, 'k');
hline = refline(0, 0);
title('GDP'); 

subplot(3,2,2); 
plot(k_uu, 'k');
hline = refline(0, 0);
title('Capital'); 

subplot(3,2,3); 
plot(hm_uu, 'k');
hline = refline(0, 0);
title('Male Market Labor Supply'); 

subplot(3,2,4); 
plot(hf_uu, 'k'); 
hline = refline(0, 0);
title('Female Market Labor Supply');


subplot(3,2,5); 
plot(lm_uu, 'k');
hline = refline(0, 0);
title('Male Home Labor Supply'); 

subplot(3,2,6); 
plot(lf_uu, 'k'); 
hline = refline(0, 0);
title('Female home Labor Supply');



%% Distributional Characteristics of the Labor: 
scatterhist(hf,hm);figure(gcf); 
xlabel('Female Market labor Supply'); 
ylabel('Male Market labor Supply'); 
%close all; 

%% Wage gap IRF:
figure; 
subplot(1,2,1); 
plot(w_gap_ee,'k'); 
title('Wage Gap Reaction to A_t Shock')
hline = refline(0, 0);
subplot(1,2,2); 
plot(w_gap_uu,'k'); 
title('Wage Gap Reaction to D_t Shock')
hline = refline(0, 0);

%close all; 
figure; 
histfit(w_gap,100,'kernel'); 
legend('Wage Gap', 'Kernel Density'); 





 