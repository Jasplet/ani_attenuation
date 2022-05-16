function [difrs, dtstars]=msac_measure_dtstar_iter(tr1,tr2,wbeg,wend,fref)
% Measure the differential t* between two traces (tr1-tr2), based on 
% equalising their windowed) weighted instantaneous frequency.

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
dtstars(2)  = 0.5; % t* attref in Mathenay and Nowack (1995)
i = 2;
while (abs(difr) >= 1e-3) && (i <=20) && (dtstars(i) < 4)
   tr_attenuated = msac_apply_tstar_operator(tr_to_attenuate, fref, dtstars(i));
   ifr = msac_ifa_wwind(tr_attenuated, wbeg, wend);
   difr = ifr - ifr_obs; % Difference in IFr between the traces
   difrs(i) = difr;
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