clear all; 
%% Steady State Solver for RBC with Household Production
% Six variables {c,hm,hf,lm,lf,k}
% Five equilibrium equations + resources constraint
% m and f denotes gender
% h is hours of market work 
% l is hours of home work

% Parameters: 
beta = 0.98; 
%lambda = 0.5;  % We will choose other value: 
alpha = 0.5; 
theta = 1/3; 
gamma = 1/3; 
delta = 0.01; 
v = 1/2; 

% Variables Order : 
% c = x(1), hm = x(2), hf = x(3), lm = x(4), lf = x(5), k = x(6)

% Solutions to steady state for many values of lambda: 
lambda = 0.1:0.01:0.9; % lambda belongs to the range 0.1-0.9:

for i = 1:1:length(lambda) ; 
    lambd(i) =  lambda(1,i); 

    F = @(x)[ 
    1-beta*(1+theta*x(6)^(theta-1)*x(2)^(gamma)*x(3)^(1-theta-gamma)-delta) ; 
    alpha*v/x(4)-(1-alpha)*lambd(i)/(1-x(2)-x(4))-alpha/x(1)*v*x(4)^(v-1)*x(5)^(1-v) ; 
    alpha*(1-v)/x(5)-(1-alpha)*(1-lambd(i))/(1-x(3)-x(5))-alpha/x(1)*(1-v)*x(4)^v*x(5)^(-v) ;
    (1-alpha)*lambd(i)/(1-x(2)-x(4))-alpha/x(1)*gamma*x(6)^(theta)*x(2)^(gamma-1)*x(3)^(1-theta-gamma); 
    (1-alpha)*(1-lambd(i))/(1-x(3)- x(5))-alpha/x(1)*(1-theta-gamma)*x(6)^(theta)*x(2)^(gamma)*x(3)^(-theta-gamma); 
    
    
   ]  
%Initial conditions: 
c0 = 1.1; 
hm0 = 0.3; 
hf0 = 0.3; 
lm0 = 0.3; 
lf0 = 0.3; 
k0 = 13.3492; % Dynare Equilibrium
% Options for the solver: 
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxFunEvals',1500);
ss = fsolve(F, [c0,hm0,hf0,lm0,lf0,k0],options);

% Additional variables:
cd = ss(1,4)^(v)*ss(1,5)^(1-v); 
U_m = alpha*(log(ss(1,1))+log(cd))+(1-alpha)*log(1-ss(1,2)-ss(1,4)); 
U_f = alpha*(log(ss(1,1))+log(cd))+(1-alpha)*log(1-ss(1,3)-ss(1,5)); 
U_hh = lambd(i)*U_m+(1-lambd(i))*U_f; 
M_ss(i,:) = [ss, cd, U_m, U_f, U_hh] ; 
end 


%% Equilibrium for each lambda: 
% This plots produces Figure 1 from the paper + others: 
subplot(3,2,1); 
plot(lambda,M_ss(:,1)); title('Steady State Comnsumption'); 
xlabel('\lambda')

subplot(3,2,2); 
plot(lambda,M_ss(:,6)); title('Steady State Capital'); 
xlabel('\lambda'); ylim([13 13.5]); 

subplot(3,2,3); 
plot(lambda,M_ss(:,2)); title('Market Labor Supply'); hold on; 
plot(lambda,M_ss(:,3)); legend('male','female'); 
xlabel('\lambda'); 

subplot(3,2,4); 
plot(lambda,M_ss(:,4)); title('Home Labor Supply'); hold on; 
plot(lambda,M_ss(:,5)); legend('male','female'); 
xlabel('\lambda'); 


subplot(3,2,5); 
plot(lambda,M_ss(:,8)); title('Individual Welfare'); hold on; 
plot(lambda,M_ss(:,9)); legend('male','female'); 
xlabel('\lambda'); 


subplot(3,2,6); 
plot(lambda, M_ss(:,10)); title('Household Welfare'); 
xlabel('\lambda'); 

%% Distributional account of shares and lambda: 
% Lambda and Home-Production: 
% Now we want to change lambda (0.1,0.9) and v (0.25,0.75)

vv = 0.3:0.01:0.6; 
for j = 1:1:length(vv); 
    v(j) = vv(1,j); 
for i = 1:1:length(lambda) ; 
    
    lambd(i) =  lambda(1,i); 
    F = @(x)[ 
    1-beta*(1+theta*x(6)^(theta-1)*x(2)^(gamma)*x(3)^(1-theta-gamma)-delta) ; 
    alpha*v(j)/x(4)-(1-alpha)*lambd(i)/(1-x(2)-x(4))-alpha/x(1)*v(j)*x(4)^(v(j)-1)*x(5)^(1-v(j)) ; 
    alpha*(1-v(j))/x(5)-(1-alpha)*(1-lambd(i))/(1-x(3)-x(5))-alpha/x(1)*(1-v(j))*x(4)^v(j)*x(5)^(-v(j)) ;
    (1-alpha)*lambd(i)/(1-x(2)-x(4))-alpha/x(1)*gamma*x(6)^(theta)*x(2)^(gamma-1)*x(3)^(1-theta-gamma); 
    (1-alpha)*(1-lambd(i))/(1-x(3)- x(5))-alpha/x(1)*(1-theta-gamma)*x(6)^(theta)*x(2)^(gamma)*x(3)^(-theta-gamma); 
    ]  ; 
c0 = 1.1; 
hm0 = 0.3; 
hf0 = 0.3; 
lm0 = 0.3; 
lf0 = 0.3; 
k0 = 13.3492; % Dynare Equilibrium
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxFunEvals',1500);
ss = fsolve(F, [c0,hm0,hf0,lm0,lf0,k0],options);
%bar(steady_state); 

% Additional vars:
cd = ss(1,4)^(v(j))*ss(1,5)^(1-v(j)); 
U_m = alpha*(log(ss(1,1))+log(cd))+(1-alpha)*log(1-ss(1,2)-ss(1,4)); 
U_f = alpha*(log(ss(1,1))+log(cd))+(1-alpha)*log(1-ss(1,3)-ss(1,5)); 
U_hh = lambd(i)*U_m+(1-lambd(i))*U_f ; 
M_ss(i,:) = [ss, cd, U_m, U_f, U_hh] ; 

LF = M_ss(:,5); 
end 

Distribution_LF(:,j) = LF; 


end 

matlab.io.saveVariablesToScript('Vars.mat', 'Distribution_LF');  
%}



% Lambda and Market production: 
ggamma = 0.1:0.01:0.5; 

for j = 1:1:length(ggamma); 
    gamm(j) = ggamma(1,j); 
for i = 1:1:length(lambda) ; 
    
    lambd(i) =  lambda(1,i); 
    F = @(x)[ 
    1-beta*(1+theta*x(6)^(theta-1)*x(2)^(gamm(j))*x(3)^(1-theta-gamm(j))-delta) ; 
    alpha*v/x(4)-(1-alpha)*lambd(i)/(1-x(2)-x(4))-alpha/x(1)*v*x(4)^(v-1)*x(5)^(1-v) ; 
    alpha*(1-v)/x(5)-(1-alpha)*(1-lambd(i))/(1-x(3)-x(5))-alpha/x(1)*(1-v)*x(4)^v*x(5)^(-v) ;
    (1-alpha)*lambd(i)/(1-x(2)-x(4))-alpha/x(1)*gamm(j)*x(6)^(theta)*x(2)^(gamm(j)-1)*x(3)^(1-theta-gamm(j)); 
    (1-alpha)*(1-lambd(i))/(1-x(3)- x(5))-alpha/x(1)*(1-theta-gamm(j))*x(6)^(theta)*x(2)^(gamm(j))*x(3)^(-theta-gamm(j)); 
    ]  ; 
c0 = 1.1; 
hm0 = 0.3; 
hf0 = 0.3; 
lm0 = 0.3; 
lf0 = 0.3; 
k0 = 13.3492; % Dynare Equilibrium
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxFunEvals',1500);
ss = fsolve(F, [c0,hm0,hf0,lm0,lf0,k0],options);
%bar(steady_state); 

% Additional vars:
cd = ss(1,4)^(v)*ss(1,5)^(1-v); 
U_m = alpha*(log(ss(1,1))+log(cd))+(1-alpha)*log(1-ss(1,2)-ss(1,4)); 
U_f = alpha*(log(ss(1,1))+log(cd))+(1-alpha)*log(1-ss(1,3)-ss(1,5)); 
U_hh = lambd(i)*U_m+(1-lambd(i))*U_f ; 
M_ss(i,:) = [ss, cd, U_m, U_f, U_hh] ; 

HF = M_ss(:,3); 
end 

Distribution_HF(:,j) = HF; 


end 



matlab.io.saveVariablesToScript('Vars2.mat', 'Distribution_HF');  



% Once variables are storaged in Vars.mat and vars2.mat, we run the code
% EXPERIMENTS_SS_RBC_Household.m; 








