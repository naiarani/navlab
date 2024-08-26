clear all;
close all;

tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];

% % Ranges of changing inputs
% deldop_n = 5:2.5:12.5; tau_n = 0.5:0.25:1; maxdop_n = 50:12.5:75;
% dopp_offset_n = 2.5:2.5:7.5; chip_delay_n = 0.05:0.05:0.15; phase_n = 30:15:90;


%% Dopp_Resolution
inputs;
for c = 35:10:55
    for deldop = 5:2.5:12.5
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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_d_o_p: %d', -5 + (i-1)*5));
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


%% Tau
inputs;
for c = 35:10:55
    for tau = 0.5:0.25:1
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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_t_a_u: %d', 0.5 + (i-1)*0.25));
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_t_a_u: %d', 0.5 + (i-1)*0.25));
        xlabel('CN0dBHz');
        ylabel('Decorrelated Variance');
        title(sprintf('SigDiag2 Column %d', col));
        grid on;
        legend('show');
    end
end

%%


%% Max_Dop
inputs;
for c = 35:10:55
    for maxdop = 50:12.5:75
        sim_ccaf_8_2;
        get_ccaf_cov_11_1;
        test_estimator_6;
        nc = [nc;xhat1];
        dc = [dc;xhat2];
        ncv = [ncv;sigdiag1];
        dcv = [dcv;sigdiag2];
        % zc = [zc;dz2];
        save(['path/maxdop_' num2str(maxdop) '.mat'],'maxdop')
    end
end

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_maxdopp_: %d', 50 + (i-1)*12.5));
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_maxdopp_: %d', 50 + (i-1)*12.5));
        xlabel('CN0dBHz');
        ylabel('Decorrelated Variance');
        title(sprintf('SigDiag2 Column %d', col));
        grid on;
        legend('show');
    end
end

%%

%% Dopp_Offset
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
        save(['path/dopp_offset' num2str(dopp_offset) '.mat'],'dopp_offset')
    end
end

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_doppoffset_: %d', 2.5 + (i-1)*2.5));
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_dopoffset_: %d', 2.5 + (i-1)*2.5));
        xlabel('CN0dBHz');
        ylabel('Decorrelated Variance');
        title(sprintf('SigDiag2 Column %d', col));
        grid on;
        legend('show');
    end
end

%%

%% Chip_Delay
inputs;
for c = 35:10:55
    for chip_delay = 0.05:0.05:0.15
        sim_ccaf_8_2;
        get_ccaf_cov_11_1;
        test_estimator_6;
        nc = [nc;xhat1];
        dc = [dc;xhat2];
        ncv = [ncv;sigdiag1];
        dcv = [dcv;sigdiag2];
        % zc = [zc;dz2];
        save(['path/chip_delay' num2str(chip_delay) '.mat'],'chip_delay')
    end
end

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_chipdelay_: %d', 0.05 + (i-1)*0.05));
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_chipdelay_: %d', 0.05 + (i-1)*0.05));
        xlabel('CN0dBHz');
        ylabel('Decorrelated Variance');
        title(sprintf('SigDiag2 Column %d', col));
        grid on;
        legend('show');
    end
end

%%


%% Phase
inputs;
for c = 35:10:55
    for phase = 30:15:90
        sim_ccaf_8_2;
        get_ccaf_cov_11_1;
        test_estimator_6;
        nc = [nc;xhat1];
        dc = [dc;xhat2];
        ncv = [ncv;sigdiag1];
        dcv = [dcv;sigdiag2];
        % zc = [zc;dz2];
        save(['path/phase' num2str(chip_delay) '.mat'],'phase')
    end
end

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_phase_: %d', 30 + (i-1)*15));
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_phase_: %d', 30 + (i-1)*15));
        xlabel('CN0dBHz');
        ylabel('Decorrelated Variance');
        title(sprintf('SigDiag2 Column %d', col));
        grid on;
        legend('show');
    end
end

%%

toc