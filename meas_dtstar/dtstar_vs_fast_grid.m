function dtstar_vs_fast_grid()

% set up parameters
noise = 0.001 ;

fast_true = 30;
tlag_true = 1.5 ;
tstar = 0.5 ;
fref = 1 ;

spol = [10, 70, 150];

fast = -90:1:90 ;
dtstar=[0:0.02:1] ;

% make the grids
[DSTARG,FASTG] = meshgrid(dtstar,fast) ;

difr=FASTG.*0 ;
ifrF=FASTG.*0 ;
ifrS=FASTG.*0 ;

% loop over SPOL
for k=1:length(spol)

   % generate synthetics
   [trN,trE,trZ]=msac_splitwave(fast_true,tlag_true,'spol',spol(k),...
      'noise',noise,'wavelet','dgauss') ;

   % apply the tstar value
   [trF,trS]=msac_rotate(trN,trE,fast_true) ;
   trSA = msac_apply_tstar_operator(trS,fref,tstar) ;

   [trN,trE]=msac_rotate(trF,trSA,-fast_true) ;

for i=1:length(fast)
   for j=1:length(dtstar)
      [trF,trS]=msac_rotate(trN,trE,fast(i)) ;
      trF = msac_apply_tstar_operator(trF,fref,dtstar(j)) ;
      
      ifrF(i,j)=msac_ifa_wwind(trF,trF.a,trF.f) ;
      ifrS(i,j)=msac_ifa_wwind(trS,trS.a,trS.f) ;
      difr(i,j)=difr(i,j)+ abs(ifrF(i,j) - ifrS(i,j)) ;
   end

end



figure
pcolor(DSTARG,FASTG,difr)
colormap(jet)
shading interp
colorbar
hold on
%contour(DSTARG,FASTG,difr,[0 0],'k-')
xlabel('dtstar')
ylabel('ref. frame rotation')
%title(sprintf('SPOL = %3.0f',spol(k)));
end



