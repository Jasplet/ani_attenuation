% MSAC_IFA_WWIND
% Calculate the (amplitude-weighted) average instantaneous frequency in a window.
%
% Usage: 
%     [infreqw] = msac_ifa(tr,wbeg,wend)                    
%         Returns the (average) instantaneous frequency of an input window,
%         weighted by the amplitude
%
% REF: Matheney, M. P. & Nowack, R. L. Geophys J Int 123, 1–15 (1995).
%

% Copyright (c) 2017-2022 James Wookey
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

function [infreqw]=msac_ifa_wwind(tr,wbeg,wend)

%  generate time base
   time = tr.b+[0:tr.npts-1]*tr.delta ;

%  locate the appropriate window
   ind = find(time>=wbeg & time<=wend) ;
   
%  generate (windowed) complex trace 
   Y = hilbert(tr.x1(ind)) ; 

%  calculate instantaneous amplitude
   inamp = sqrt( real(Y).^2 + imag(Y).^2 ) ;

%  instantaneous phase
   inph = atan(imag(Y)./real(Y)) ;

%  calculate instantaneous frequency trace (Matheney equ. 7)
   yy = real(Y);
   ys = imag(Y);

   dyydt = gradient(yy,tr.delta) ;
   dysdt = gradient(ys,tr.delta) ;

   infreq = (1/(2*pi)) .* (yy.*dysdt-ys.*dyydt)./(yy.^2+ys.^2) ;
   
%  calculate amplitude weighted average (Matheney equ. 9)
   infreqw = sum(infreq.*inamp.^2) ./ sum(inamp.^2) ;

end