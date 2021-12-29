%% Mod File: RBC With Household Production: 

% VARIABLES: 
var c cd hm lm hf lf k A z D eta w_gap H Y ;

% EXOGENOUS VARIABLES: 
varexo ee uu; 

% PARAMETERS: 
parameters alppha rrho thetta betta lammbda ggamma deltta v xi; 

% CALIBRATION: 
betta = 0.98; 
rrho = 0.9;
xi = 0.9;  
deltta = 0.01; 
thetta = 1/3; 
ggamma = 1/3; 
alppha = 1/2; 
lammbda = 0.6; 
v = 0.5; 


% EQUILIBRIUM: 

model; 

% productivity: 
z = rrho*z(-1)+ee;
eta = xi*eta(-1)+uu; 
A = exp(z); 
D = exp(eta); 

% resource constraint: 
c + D*lm^v*lf^(1-v) + k = A*k(-1)^(thetta)*hm^(ggamma)*hf^(1-thetta-ggamma)+(1-deltta)*k(-1); 

% household decisions
% Eurler eq: 
(1/c) = betta*(1+A*thetta*k^(thetta-1)*(hm(+1))^(ggamma)*(hf(+1))^(1-thetta-ggamma)-deltta)*(1/c(+1)) ; 

% (lm): 
alppha*v / lm = (1-alppha)*lammbda/(1-hm-lm)+(alppha/c)*D*v*lm^(v-1)*lf^(1-v); 

% (lf): 
alppha*(1-v) / lf = (1-alppha)*(1-lammbda)/(1-hf-lf)+(alppha/c)*D*(1-v)*lm^v*lf^(-v); 

% (hm): 
(1-alppha)*lammbda / (1-hm-lm) = (alppha/c)*ggamma*A*k(-1)^(thetta)*hm^(ggamma-1)*hf^(1-ggamma-thetta); 

% (hf): 
(1-alppha)*(1-lammbda) / (1-hf-lf) = (alppha/c)*(1-thetta-ggamma)*A*k(-1)^(thetta)*hm^(ggamma)*hf^(-ggamma-thetta); 

% (cd)
cd = lm+lf; 

w_gap = ggamma / (1-thetta-ggamma) * hf/hm;

H = hm + hf; 
Y = A*k(-1)^(thetta)*hm^(ggamma)*hf^(1-thetta-ggamma)+(1-deltta)*k(-1); 
end; 


% INITIAL VALUES: 
% Computed in: Steady_State_RBC_household_production.m
initval; 
c  = 0.9475; 
hm = 0.27; 
hf = 0.27; 
lm = 0.24; 
lf = 0.24; 
k  = 10; 
cd = lm+lf; 
z = 0; 
A = 1;
eta = 0;  
D = 1 ;
w_gap = 1; 
end; 

% PRINT STEADY STATE: 
steady; 
check; 

% SHOCKS: 
shocks; 
var ee; stderr 0.02; 
var uu; stderr 0.02;
end;

%/////////////////////////////////////////////
% SOLVING THE MODEL: 

% STOCHASTIC SIMULATION: 

stoch_simul(order=2,irf=100,dr_algo=0,periods=1000); 


% DETERMINISTIC SIMULATION: 
%perfect_foresight_setup(periods=400);
%perfect_foresight_solver(stack_solve_algo=0);





















 