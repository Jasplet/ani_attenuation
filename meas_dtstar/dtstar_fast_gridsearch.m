function [difr] = dtstar_fast_gridsearch(trN, trE, wbeg, wend, s_fast, s_dtstar, fref)
% Performs a grid search over fast direction and differential attenuation 
n = length(s_fast);
m = length(s_dtstar);
ifrF = zeros(n,m);
ifrS = zeros(n,m);
difr = zeros(n,m); 
for i=1:length(s_fast)
   for j=1:length(s_dtstar)
      [trF,trS]=msac_rotate(trN,trE,s_fast(i)) ;
      trF = msac_apply_tstar_operator(trF,fref,s_dtstar(j)) ;
      
      ifrF(i,j)=msac_ifa_wwind(trF,wbeg,wend) ;
      ifrS(i,j)=msac_ifa_wwind(trS,wbeg,wend) ;
      difr(i,j)=abs(ifrF(i,j) - ifrS(i,j)) ;
   end

end
end

