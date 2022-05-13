% measuring dt* test.
clc; clear ; close all
% set up parameters
noise = 0.001 ;
spol = -75;
fast_true = 40;
tlag_true = 0 ;
tstar = 1;
fref=10 ;

% generate synthetics
[trN,trE,trZ]=msac_splitwave(fast_true,tlag_true,spol,noise) ;
samps = [1:length(trN.x1)]*trN.delta;

% apply the tstar value
[trF,trS]=msac_rotate(trN,trE,fast_true) ;
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;

% Try to measure t*
ifrtrF=msac_ifa_wwind(trF,trF.a,trF.f) ;
ifrtrS=msac_ifa_wwind(trSA,trSA.a,trSA.f) ;

[difrs, dtstars] =  meas_tstar(trF, trSA, trF.a, trF.f, fref);

srange = [0:0.2:4];
difr_range = zeros(1,length(srange));
for i = 1:length(srange)
   tr_attenuated = msac_apply_tstar_operator(trF, fref, srange(i));
   ifr_trial = msac_ifa_wwind(tr_attenuated, tr_attenuated.a, tr_attenuated.f);
   difr = ifr_trial - ifrtrS;
   difr_range(i) = difr;
end

% Plot iterative search for t* and search range
figure() 
plot(srange, difr_range, 'k-', 'Linewidth',1.5)
hold on
plot(dtstars, difrs, 'rx', 'Markersize',10)
xlabel('\Delta t*')
ylabel('\Delta IFr')

%% Test measurement method on a real world example from YKW3

trE = msac_read('data/YKW3_splitex.BHE');
trN = msac_read('data/YKW3_splitex.BHN');
trZ = msac_read('data/YKW3_splitex.BHZ');
wbeg = trE.a;
wend = trE.f;

fast_true = 44;
spol = rem(148, 180) - 90 ;
tlag_true = 0.95;
fref = 10;

[trF,trS]=msac_rotate(trN,trE,fast_true) ;

ifrF=msac_ifa_wwind(trF, wbeg, wend) ;
ifrS=msac_ifa_wwind(trS, wbeg, wend) ;

difr = ifrF - ifrS;

[difrs, dtstars] =  meas_tstar(trF, trS, wbeg, wend, fref);

srange = [0:0.2:4];
difr_range = zeros(1,length(srange));
for i = 1:length(srange)
   tr_attenuated = msac_apply_tstar_operator(trF, fref, srange(i));
   ifr_trial = msac_ifa_wwind(tr_attenuated, tr_attenuated.a, tr_attenuated.f);
   difr = ifr_trial - ifrS;
   difr_range(i) = difr;
end

% Plot iterative search for t* and search range
figure() 
plot(srange, difr_range, 'k-', 'Linewidth',1.5)
hold on
plot(dtstars, difrs, 'rx', 'Markersize',10)
xlabel('\Delta t*')
ylabel('\Delta IFr')
%%
function [difrs, dtstars] =  meas_tstar(tr1, tr2, wbeg, wend, fref)

ifrtr1=msac_ifa_wwind(tr1, wbeg, wend) ;
ifrtr2=msac_ifa_wwind(tr2, wbeg, wend) ;

difr = ifrtr1 - ifrtr2 ;

% based on sign, determine direction of correction
signTS = sign(difr) ;
if signTS<0
   tr_to_attenuate = tr2 ; % This is (confisuingly) the "reference pulse" in Matheny and Nowack (1995)
   ifr_obs = ifrtr1 ; % this is f_obs in Mathenay and Nowack (1995)
   ifr_old = ifrtr2 ; % f_attn in Matheney and Nowack for initial condition t* = 0 
else
   tr_to_attenuate = tr1 ;
   ifr_obs = ifrtr2 ;
   ifr_old = ifrtr1 ; % f_attn in Matheney and Nowack for initial condition t* = 0 
end
difrs = zeros(50,1);
dtstars = zeros(50,1);
% Initial conditions. t* = 0
dtstars(1) = 0;
difrs(1) = difr;
% Start search at t* = 0.05
dtstars(2)  = 0.05; % t* attref in Mathenay and Nowack (1995)
i = 2;
while (abs(difr) >= 1e-3) && (i <=20)
   tr_attenuated = msac_apply_tstar_operator(tr_to_attenuate, fref, dtstars(i));
   ifr = msac_ifa_wwind(tr_attenuated, wbeg, wend);
   difrs(i) = ifr - ifr_obs; % Difference in IFr between the traces
   step = dtstars(i) - dtstars(i-1);
   dfdts = (ifr - ifr_old)/step;
   ifr_old = ifr;
   % set up next iteration
   dtstars(i+1) = dtstars(i) - difrs(i)/dfdts;
   i = i+1;
   
end
difrs = difrs(1:i);
dtstars = dtstars(1:i).*signTS;
end


