clear all;
close all;

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];
colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};


% % Ranges of changing inputs
% deldop_n = 5:2.5:12.5; tau_n = 0.5:0.25:1; maxdop_n = 50:12.5:75;
% dopp_offset_n = 2.5:2.5:7.5; chip_delay_n = 0.05:0.05:0.15; phase_n = 30:15:90;


%% Max Doppler

tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];


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
    end
end

toc
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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_m_a_x_d_o_p: %.2f', 50 + (i-1)*12.5));
        xlabel('CN0-dBHz');
        ylabel('Correlated $\sigma$', 'Interpreter','latex');
        title(sprintf('Correlated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_m_a_x_d_o_p: %.2f', 50 + (i-1)*12.5));
        xlabel('CN0-dBHz');
        ylabel('Decorrelated $\sigma$', 'Interpreter','latex');
        title(sprintf('Decorrelated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
    end
end

%%


%% Tau

tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];
colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};

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
toc
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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_m_a_x_t_a_u: %.2f', 0.5 + (i-1)*0.25));
        xlabel('CN0-dBHz');
        ylabel('Correlated $\sigma$', 'Interpreter','latex');
        title(sprintf('Correlated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');;
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_m_a_x_t_a_u: %.2f', 0.5 + (i-1)*0.25));
        xlabel('CN0-dBHz');
        ylabel('Decorrelated $\sigma$', 'Interpreter','latex');
        title(sprintf('Decorrelated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
    end
end

%%


%% Dop_Resolution

tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];
colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};


inputs;
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
toc

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_d_o_p_r_e_s: %.2f', 5 + (i-1)*2.5));
        xlabel('CN0-dBHz');
        ylabel('Correlated $\sigma$', 'Interpreter','latex');
        title(sprintf('Correlated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_d_o_p_r_e_s: %.2f', 5 + (i-1)*2.5));
        xlabel('CN0-dBHz');
        ylabel('Decorrelated $\sigma$', 'Interpreter','latex');
        title(sprintf('Decorrelated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
    end
end

%%

%% Dopp_Offset
tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];
colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};


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

toc

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_d_o_p_o_f_f_s_e_t_: %.2f', 2.5 + (i-1)*2.5));
        xlabel('CN0-dBHz');
        ylabel('Correlated $\sigma$', 'Interpreter','latex');
        title(sprintf('Correlated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_d_o_p_o_f_f_s_e_t_: %.2f', 2.5 + (i-1)*2.5));
        xlabel('CN0-dBHz');
        ylabel('Decorrelated $\sigma$', 'Interpreter','latex');
        title(sprintf('Decorrelated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
    end
end

%%

%% Chip_Delay
tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];
colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};


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
    end
end
toc

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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_c_h_i_p_d_e_l_a_y_: %.2f', 0.05 + (i-1)*0.05));
        xlabel('CN0-dBHz');
        ylabel('Correlated $\sigma$', 'Interpreter','latex');
        title(sprintf('Correlated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_c_h_i_p_d_e_l_a_y_: %.2f', 0.05 + (i-1)*0.05));
        xlabel('CN0-dBHz');
        ylabel('Decorrelated $\sigma$', 'Interpreter','latex');
        title(sprintf('Decorrelated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
    end
end

%%


%% Phase
tic

% add spoof chip center, dopp center, phase

nc = [];
dc = [];
ncv = [];
dcv = [];
zc = [];
colnames = {'True Amplitude', 'True Tau', 'True Doppler', 'True Phase', 'Spoofed Amplitude', 'Spoofed Tau', 'Spoofed Doppler', 'Spoofed Phase'};

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
    end
end
toc
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
        plot(35:10:55, submatrices_ncv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_p_h_a_s_e_: %.0f', 30 + (i-1)*15));
        xlabel('CN0-dBHz');
        ylabel('Correlated $\sigma$', 'Interpreter','latex');
        title(sprintf('Correlated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
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
        plot(35:10:55, submatrices_dcv{i}(:, col), 'o-', 'DisplayName', sprintf('ncv_p_h_a_s_e_: %.0f', 30 + (i-1)*15));
        xlabel('CN0-dBHz');
        ylabel('Decorrelated $\sigma$', 'Interpreter','latex');
        title(sprintf('Decorrelated %s', colnames{col}));
        grid on;
        legend('Location', 'northoutside');
    end
end

%%
