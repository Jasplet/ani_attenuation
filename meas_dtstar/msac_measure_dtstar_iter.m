function [difrs, dtstars]=msac_measure_dtstar_iter(tr1,tr2,wbeg,wend, dtstar_init, fref)
% Measure the differential t* between two traces (tr1-tr2), based on 
% equalising their windowed) weighted instantaneous frequency.

ifrtr1=msac_ifa_wwind(tr1, wbeg, wend) ;
ifrtr2=msac_ifa_wwind(tr2, wbeg, wend) ;
difr_test = ifrtr1 - ifrtr2 ;

% based on sign, determine direction of correction
signTS = sign(difr_test) ;
if signTS<0
   tr_to_attenuate = tr2 ; % This is (confisuingly) the "reference pulse" in Matheny and Nowack (1995)
   ifr_obs = ifrtr1 ; % this is f_obs in Mathenay and Nowack (1995)
   ifr_ref = ifrtr2 ; % f_attn in Matheney and Nowack for initial condition t* = 0 
else
   tr_to_attenuate = tr1 ;
   ifr_obs = ifrtr2 ;
   ifr_ref = ifrtr1 ; % f_attn in Matheney and Nowack for initial condition t* = 0 
end

% Gradient descent method
% Initial conditions. t* = 0
ifr_old = ifr_ref;
difrs = zeros(50,1);
dtstars = zeros(50,1);
difr = ifr_ref - ifr_obs;
difrs(1) = difr;
% Start search at initial guess t* (this gives us a second points along
% with t*  = 0)
dtstars(2)  = dtstar_init; % t* attref in Mathenay and Nowack (1995)
i = 2;
while (abs(difr) >= 1e-3)
   tr_attenuated = msac_apply_tstar_operator(tr_to_attenuate, fref, dtstars(i));
   ifr = msac_ifa_wwind(tr_attenuated, wbeg, wend);
   difr = ifr - ifr_obs; % Difference in IFr between the traces
   difrs(i) = difr;
   step = dtstars(i) - dtstars(i-1);
   dfdts = (ifr - ifr_old)/step;
   ifr_old = ifr;
   % set up next iteration
   dtstars(i+1) = dtstars(i) - difr/dfdts;
   i = i+1;
   % Try to catch any cases wherer iteration breaks down (e.g. bad
   % waveforms where difr != 0)
   if (dtstars(i+1) > 4.0) || (i > 20)
       fprintf('Warning bad dt* measurement. dIFr not approaching 0 within a reasonable search range\n')
       %dtstars(i+1) = NaN;
       break
   end
end
difrs = difrs(1:i);
dtstars = dtstars(1:i).*signTS;

end
