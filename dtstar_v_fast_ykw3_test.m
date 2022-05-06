function [ ]= dtstar_v_fast_ykw3_test()
% Testing a revised version of James' dtstar v fast function but using a
% real SKS result rather than a synthetic

% set up parameters
fast_meas = 58.0;
spol = rem(108.725, 180) - 90 ;
tlag_meas = 1.10;
fref = 25;
tstar = 0.05;


trE = msac_read('ykw3_test.BHE');
trN = msac_read('ykw3_test.BHN');
trZ = msac_read('ykw3_test.BHZ');
wbeg = trE.a;
wend = trE.f;

fast = -90:90 ;

%trSA = trS ;
for i=1:length(fast)
   [trF,trS]=msac_rotate(trN,trE,fast(i)) ;
   trS = msac_tshift(trS,-tlag_meas,'int') ;
   ifrF(i)=msac_ifa_wwind(trF,wbeg,wend) ;
   ifrS(i)=msac_ifa_wwind(trS,wbeg,wend) ;
   dtstar(i) = msac_measure_dtstar(trF,trS,wbeg,wend,fref,[0:0.02:5]) ;
   
   M = cov(trF.x1,trS.x1) ;
   E = eig(M) ;

   lam1(i) = max(E);    
   lam2(i) = min(E) ;
end
%% Plot dt* v fast direction
figure("Position",[0 0 600 2000]) ;
plot(fast,dtstar,'b-','LineWidth',1.5) ;
hold on 
aax=axis() ;

plot([fast_meas fast_meas],[-100 100],'k-','LineWidth',1.5) ;
plot([spol spol],[-100 100],'k--','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('dt* (s)')
title(['tlag = ', sprintf('%5.2fs',tlag_meas), '; dt* =', sprintf('%5.2fs',tstar)])

axis([fast(1) fast(end) -5 5]) ;

%% Plot dIfr v fast
figure()
difr = ifrF - ifrS;
plot(fast,difr,'r-','LineWidth',1.5) ;
xline(fast_meas,'-k','LineWidth',1.5);
xline(spol, '--k','LineWidth', 1.5);
aax=axis() ;
xlim([tast(1), fast(end)]);
xlabel('Reference frame rotation (degrees)')
ylabel('Difference in Inst. freq. (Hz)')
legend('difr')
title(['tlag = ' sprintf('%5.2fs',tlag_meas) '; dt* =' sprintf('%5.2fs',tstar)])

end

