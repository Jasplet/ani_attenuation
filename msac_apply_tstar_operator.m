% msac_apply_tstar_operator - attenuate a seismic trace using a causal
%                             t* operator.
%
% [tr_out] = msac_apply_tstar_operator(tr,fref,tstar)
%
% Usage: 
%     [tr_out] = msac_apply_tstar_operator(tr,fref,tstar)                    
%         Apply a causal t* operator (tstar, at a reference frequency fref)
%         to trace tr, returning its attenuated equivalent (tr_out)
%
% REF: Matheney, M. P. & Nowack, R. L. Geophys J Int 123, 1–15 (1995).
%
% 2022 version with proper negative frequency handling.
%
function [tr_out] = msac_apply_tstar_operator(tr,fref,tstar)

%  extract some useful information
   t = [0:tr.npts-1].*tr.delta ;
   y = tr.x1 ;
   Fs = 1/tr.delta ;

%  zero pad trace then FFT
   n = 2^nextpow2(tr.npts) ;
   Y = fft(y,n) ;

%  frequency vector
   f = Fs*(0:(n/2))/n ;

%  angular frequency
   om = 2.*pi.*f ;
   omref = 2*pi*fref ;

%  create (causal) t* multiplier
   D = zeros(1,n/2+1) ;
   D(1) = complex(1,0) ;
   for ii=2:1+n/2
      D(ii) = exp ( complex( -(1/2)*om(ii)*tstar , ...
                              (1/pi)*om(ii)*tstar*log(om(ii)/omref) ) ) ;
   end

%  apply to FD signal
   for ii = 2:1+n/2
      Y(ii) = Y(ii).*D(ii) ;
   end
   
%  restore symmetry
   Y(2+n/2:length(Y)) = conj(Y(n/2:-1:2)) ;
      
%  inverse FFT
   X = ifft(Y) ;

%  generate output trace
   tr_out = tr ;
   tr_out.x1 = real(X(1:tr.npts)) ;
end


