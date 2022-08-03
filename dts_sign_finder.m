%% Script to test determining the sign of dt*
close all ; clear
stat = 'YIRG';
file = sprintf('/Users/ja17375/Projects/DeepMelt/Ethiopia/%s_sks_reprocessed.txt',stat);
data = readtable(file, 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = sprintf('/Users/ja17375/Projects/DeepMelt/Ethiopia/data/%s/run/', stat);

snr_min = 1;
idxs = find(data.SNR > snr_min);
n = length(idxs);
fast_meas = 16;
dts_best = 0.8;
difrs = zeros(n,1);

for i=1:n
    idx = idxs(i);
    date = char(data.DATE(idx));
    time = char(data.TIME(idx));
    wbeg = data.WBEG(idx);
    wend = data.WEND(idx);
    spol = data.SPOL(idx);
    fileid = sprintf('%s_%s_%s',stat, date, time(1:4));
    fprintf([path, fileid, '\n'])
    trN = msac_read([path, fileid, '.BHN']);
    trE = msac_read([path, fileid, '.BHE']);
    trZ = msac_read([path, fileid, '.BHZ']);
    pts = 0:trN.npts-1;
    times = trN.b + pts*trN.delta;
    ind = find((times >=wbeg) & (times <=wend));
    [trF, trS] = msac_rotate(trN, trE, fast_meas);
    ifrF = msac_ifa_wwind(trF,wbeg, wend) ;
    ifrS = msac_ifa_wwind(trS,wbeg, wend) ; 
    difrs(i) = ifrF - ifrS;
    figure()
    plot(times(ind),trF.x1(ind), 'b-')
    hold on 
    plot(times(ind), trS.x1(ind), 'r-')

    legend('"Fast"', '"Slow"')
    xlabel('Time [s]')
    title(['Station ',stat, '. Date: ',date, '-', time])
    saveas(gcf,['/Users/ja17375/Projects/DeepMelt/Ethiopia/Figs/',fileid,'_waveform.png'])
    hold off
end