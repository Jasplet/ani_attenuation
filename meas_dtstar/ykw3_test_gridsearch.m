%% Script to test grid search measurments on YKW3 data.
close all ; clear
ykw3_data = readtable('YKW3_sks_splitting.txt', 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = '/Users/ja17375/Projects/DeepMelt/CanadianShield/YKW3/run/';
% find data that has a high SNR (for now, maybe this can be relaxed in the
% future)
idx_test = find(ykw3_data.SNR > 20);
fref=10;
fasts = -90:90;
dts = 0:0.05:4;
% make the grids
[DSTARG,FASTG] = meshgrid(dts,fasts) ;
difr = FASTG.*0 ;
%idx_test = [67, 94,167,186];
difrs = zeros(length(fasts),length(dts), length(idx_test));
for i=1:length(idx_test)
    ind = idx_test(i);
    date = char(ykw3_data.DATE(ind));
    time = char(ykw3_data.TIME(ind));
    wbeg = ykw3_data.WBEG(ind);
    wend = ykw3_data.WEND(ind);
    spol = ykw3_data.SPOL(ind);

    fileid = sprintf('YKW3_%s_%s',date, time);

    trN = msac_read([path, fileid, '.BHN']);
    trN.a = wbeg;
    trN.f = wend;
    trE = msac_read([path, fileid, '.BHE']);
    trE.a = wbeg;
    trE.f = wend;
    trZ = msac_read([path, fileid, '.BHZ']);
  
    difrs(:,:,i) = dtstar_fast_gridsearch(trN, trE, fasts, dts, fref);
    difr = difr + difrs(:,:,i);
    % Plot individual measurement surface
%     figure()
%     pcolor(DSTARG,FASTG,difrs(:,:,i))
%     colormap(jet)
%     shading interp
%     colorbar
%     hold on
%     %contour(DSTARG,FASTG,difr,[0 0],'k-')
%     xlabel('dtstar')
%     ylabel('ref. frame rotation')
%     title(sprintf('YKW3. SPOL = %3.0f',spol));
end

figure
pcolor(DSTARG,FASTG,difr./length(idx_test))
colormap(jet)
shading interp
colorbar
hold on
%contour(DSTARG,FASTG,difr,[0 0],'k-')
xlabel('dtstar')
ylabel('ref. frame rotation')
title('YKW3 stacked')
