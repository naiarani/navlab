clear all;
close all;

tic

% add spoof chip center, dopp center, phase

inputs;

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];

% Ranges of changing inputs


%% Dopp_Resolution
for c = 35:10:55
    for dop_res = 5:2.5:12.5 
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
numrows_ncv = size(ncv, 1);
submatrices_ncv = cell(1, numrows_ncv/3);

for i = 1:numrows_ncv/3
    start_row = (i - 1) * 3 + 1;
    end_row = start_row + 2;     
    submatrices_ncv{i} = ncv(start_row:end_row, :);
end

num_submatrices = numel(submatrices_ncv);
figure;
for i = 1:num_submatrices
    for col=1:size(ncv,2)
        subplot(2,4,col);
        hold on;
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv: %d', -5 + (i-1)*5));
        xlabel('CN0dBHz');
        ylabel('Correlated Variance');
        title(sprintf('SigDiag1 Column %d', col));
        grid on;
        legend('show');
    end
end

numrows_dcv = size(dcv, 1);
submatrices_dcv = cell(1, numrows_dcv/3);

for i = 1:numrows_dcv/3
    start_row = (i - 1) * 3 + 1;
    end_row = start_row + 2;     
    submatrices_dcv{i} = dcv(start_row:end_row, :);
end

num_submatrices = numel(submatrices_dcv);
figure;
for i = 1:num_submatrices
    for col=1:size(dcv,2)
        subplot(2,4,col);
        hold on;
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_d_o_p: %d', -5 + (i-1)*5));
        xlabel('CN0dBHz');
        ylabel('Decorrelated Variance');
        title(sprintf('SigDiag2 Column %d', col));
        grid on;
        legend('show');
    end
end

%%

toc
