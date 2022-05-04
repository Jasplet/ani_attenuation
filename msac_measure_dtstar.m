function [dtstar]=msac_measure_dtstar(tr1,tr2,t0,t1,fref,srange)
% Measure the differential t* between two traces (tr1-tr2), based on 
% equalising their windowed) weighted instantaneous frequency.

ifrtr1=msac_ifa_wwind(tr1,t0,t1) ;
ifrtr2=msac_ifa_wwind(tr2,t0,t1) ;

difr = ifrtr1 - ifrtr2 ;

% based on sign, determine direction of correction
signTS = sign(difr) ;
if signTS<0
   tr_to_attenuate = tr2 ; % This is (confisuingly) the "reference pulse" in Matheny and Nowack (1995)
   ref_infreq = ifrtr1 ; % this is f_obs in Mathenay and Nowack (1995)
else
   tr_to_attenuate = tr1 ;
   ref_infreq = ifrtr2 ;
end

% dtstar_old = 0; % t* attref in Mathenay and Nowack (1995)
% i = 0;
% while abs(difr) > 0.1
%    tr_attenuated = msac_apply_tstar_operator(tr_to_attenuate, fref, dtstar_old);
%    f_attref = msac_ifa_wwind(tr_attenuated, t0, t1);
%    
% end

ifa = zeros(1,length(srange)) ;

% set range of tstar, and search to get ifa as a function of t*
for its=1:length(srange)
    % Attenuate reference trace by candidate value of t*
	tr_ts = msac_apply_tstar_operator(tr_to_attenuate,fref,srange(its)) ;
    % measure instantaneous frequency, ifa, of attenuated trace
    ifa(its) = msac_ifa_wwind(tr_ts,t0,t1) ;
end

% calculate tstar between fast, slow. Assign sign based on delta(IFr).
dtstar = signTS.*interp1(ifa,srange,ref_infreq,'pchip',NaN) ;

if isnan(dtstar), warning('DTSTAR out of range'), end

% fprintf('T* = %8.4f\n',dtstar) ;
% if signTS<0
%    fprintf('Negative number means that fast wave is more attenuated.\n')
% else
%    fprintf('Positive number means that slow wave is more attenuated.\n')
% end
% fprintf('\n') ;