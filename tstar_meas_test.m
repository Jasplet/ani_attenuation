% measuring dt* test.

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
[trR, trT] = msac_rotate(trN, trE, fast_true);

% apply the tstar value
[trF,trS]=msac_rotate(trN,trE,fast_true) ;
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;

% Try to measure t*
ifrtr1=msac_ifa_wwind(trF,trF.a,trF.f) ;
ifrtr2=msac_ifa_wwind(trSA,trSA.a,trSA.f) ;

difr = ifrtr1 - ifrtr2;

dtstar_old  = 0; % t* attref in Mathenay and Nowack (1995)
i = 0;
while difr < 0.05
   if i == 0
      % special case to get is going try t* = 0.05.
      dtstar = 0.05;
   else
      dtstar = dtstar_old + (difr/dfdts);
   end
   tr_attenuated = msac_apply_tstar_operator(tr1, fref, dtstar);
   ifr_trial = msac_ifa_wwind(tr_attenuated, tr_attenuated.a, tr_attenuated.f);
   difr = ref_infreq - f_attref; % Difference in IFr between the traces
   dfdtstar = difr - difr_old; 
   dtstar = dtstar_old + difr/dfdtstar;
   
end