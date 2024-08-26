
%% Correct Measurements - Dopp_Offset

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];

inputs;
for c = 35:10:55
    for dopp_offset = 2.5:2.5:7.5
        sim_ccaf_8_2;
        get_ccaf_cov_11_1;
        test_estimator_6;
        nc = [nc;xhat1];
        dc = [dc;xhat2];
        ncv = [ncv;sigdiag1];
        dcv = [dcv;sigdiag2];
        % zc = [zc;dz2];
    end
end

%%

colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};
plotData(ncv, colnames, 'Correlated $\sigma$');
plotData(dcv, colnames, 'Decorrelated $\sigma$');

toc

%%