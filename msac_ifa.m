% MSAC_IFA
% Calculate the instantaneous frequency, amplitude and phase of an input 
% seismic trace.
%
% Usage: 
%     [inamp,infreq] = msac_ifa(tr)                    
%         Produce instantaneous amplitude and frequency traces from an input
%         seismic trace (in SAC format). 
%
%	  [...] = msac_ifa(...,'weight',wlength) ;
%		  Apply an amplitude weighting to the instantaneous frequency, using
%         a sliding window of wlength (seconds). The weight is determined
%         from the envelope function.
%
%	  [inamp,infreq,inph] = msac_ifa(...)
%         Optionally also return the instantaneous phase.
%
% REF: Matheney, M. P. & Nowack, R. L. Geophys J Int 123, 1–15 (1995).
%

% Copyright (c) 2017-2021 James Wookey
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [varargout]=msac_ifa(tr,varargin)

%  weight 
   iweight = 0 ;

%  process the optional arguments
   iarg = 1 ;
   while iarg <= (length(varargin))
      switch lower(varargin{iarg})
         case 'weight'
         	iweight = 1;
            wlength = varargin{iarg+1} ;
            iarg = iarg + 2 ;
         otherwise 
            error(['Unknown option: ' varargin{iarg}]) ;   
      end   
   end 

%  first, calculate instantaneous amplitude
   Y = hilbert(tr.x1) ; 

   inamp = sqrt( real(Y).^2 + imag(Y).^2 ) ;
   infreq = 0 ;

   inph = atan(imag(Y)./real(Y)) ;

   yy = real(Y);
   ys = imag(Y);

   dyydt = gradient(yy,tr.delta) ;
   dysdt = gradient(ys,tr.delta) ;

   infreq = (1/(2*pi)) .* (yy.*dysdt-ys.*dyydt)./(yy.^2+ys.^2) ;
   infreqw = infreq.*NaN ;
   if iweight
      iwhalf = round(wlength/2/tr.delta) ;
      for i=1:length(infreq)
         i0 = max([1 i-iwhalf]) ;
         i1 = min([length(infreq) i+iwhalf]) ;
         infreqw(i) = sum(infreq(i0:i1).*inamp(i0:i1).^2) ./ sum(inamp(i0:i1).^2) ;
      end
      infreq = infreqw ;
   end

%  assign outputs
   if nargout==2
      varargout = {inamp,infreq}  ;
   elseif nargout==3
      varargout = {inamp,infreq,inph}  ;      
   else
      error('Wrong number of outputs') ;
   end
end