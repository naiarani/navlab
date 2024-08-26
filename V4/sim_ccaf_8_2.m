% clear
% close all

% rng(4);                 % random number generator seed

% INPUTS
% Simulation control parameters
Tcoh = 20;                                            %  coherent integration time (msec)
Tcohs = Tcoh/1000;                             %  coherent integration time (sec)
chip_resolution = 0.02;                           %  chips  (50 samples per chip, so 1/50 = 0.02 is smallest)
dopp_resolution = dop_res;              %  Hz  (there is minimal benefit in setting this <1/(4*Tcoh)
max_del_tau = tau;                               %  chips
max_del_doppler = maxdop; %1/Tcohs;                %  Hz
chip_center = c_center;                                %  chips
dopp_center = d_center;                            %  Hz
svs = 1;                                                %  increase to account for PRN code cross-correlation
 
% True signal
true_amplitude = 1;
true_chip_delay = chip_center + 0.0;      %  chips
true_dopp_offset = dopp_center - dopp_offset;     %  Hz
true_phase = 0*pi/180;                          % rad
CN0dBHz = c;                                      %  dB-Hz

% Spoofed signal
spoofed_amplitude = 0.9;
spoofed_chip_delay = chip_center + chip_delay;                      %  chips
spoofed_dopp_offset = dopp_center + dopp_offset;         %  Hz
spoofed_phase = phase*pi/180;                                          % rad
prn_spoofed = 17;                                                          % works for any 1-32

% multipath signal4
multipath_amplitude = 0.4*0; 
multipath_chip_delay = 0.1;                               %  chips
multipath_dopp_shift = dopp_offset;                                  %  Hz
multipath_phase = (2*rand-1)*pi;  % rad

% Load PRN codes
load CA_code_table_50MHz.mat
fs = 50e6;                  %  sample rate (Hz)
Ts = 1/fs;                   % sample interval (sec)
t = 0 : Ts : Tcoh * 1e-3 - Ts;
tmp_new_prn_order = randperm(32);
new_prn_order = [prn_spoofed tmp_new_prn_order(tmp_new_prn_order ~=prn_spoofed) ];
sv_index_spoofed = 1;
caCodesTable(1:end,1:end) = caCodesTable(new_prn_order,1:end);
codes = caCodesTable(1:svs,1:end);
codes = repmat(codes,1,Tcoh);

% Define chip spacings and Dopplers to sample 
chip_frac = chip_center + [-max_del_tau:chip_resolution:max_del_tau];
tau_corr = chip_frac'*1e-6;              % sec
i_tau_corr = round(tau_corr/Ts);      % code table column index corresponding to tau_corr
dopp_corr = dopp_center + [-max_del_doppler:dopp_resolution:max_del_doppler];       % Hz

% Randomly shift code and carrier phases for all PRNs except 'sv_index_spoofed'
rand_tau_shift = (2*rand(svs,1) - 1)*length(codes);
rand_tau_shift(sv_index_spoofed) = 0;
rts = round(rand_tau_shift);
for k = 1:svs
    codes_in(k,:) = circshift(codes(k,:),rts(k));
end
codes_in(sv_index_spoofed,:) = circshift(codes(sv_index_spoofed,:),round(true_chip_delay*1e-6/Ts));
rand_dopp_shift = (2*rand(svs,1) - 1)*max_del_doppler*10;
rand_dopp_shift(sv_index_spoofed) = 0;
carrs_in = exp(1i*(2*pi*(true_dopp_offset+rand_dopp_shift)*t + true_phase));

% Generate true signals in
true_signals_in = true_amplitude*carrs_in.*codes_in; 

% Generate spoofed signals -- assumes PRNs other than 'sv_index_spoofed' are perfectly replicated by spoofer
spoofed_carrs_in = carrs_in;
spoofed_carrs_in(sv_index_spoofed,:) = exp(1i*(2*pi*spoofed_dopp_offset*t + spoofed_phase));
spoofed_codes_in = codes_in;
spoofed_codes_in(sv_index_spoofed,:) = circshift(codes(sv_index_spoofed,:),round(spoofed_chip_delay*1e-6/Ts));
spoofed_signals_in = spoofed_amplitude*spoofed_carrs_in.*spoofed_codes_in;

% Generate multipath signal  -- only on 'sv_index_spoofed'
multipath_carrs_in = carrs_in;
multipath_carrs_in(sv_index_spoofed,:) = carrs_in(sv_index_spoofed,:).*exp(1i*(2*pi*multipath_dopp_shift*t + multipath_phase));
multipath_codes_in = codes_in;
multipath_codes_in(sv_index_spoofed,:) = circshift(codes(sv_index_spoofed,:),round(multipath_chip_delay*1e-6/Ts));
multipath_signals_in = multipath_amplitude*multipath_carrs_in.*multipath_codes_in;

% Total signal in
signal_in = true_signals_in + spoofed_signals_in + multipath_signals_in;

% Add noise
N0 = 1/(10^(CN0dBHz/10));
sig = sqrt(N0/Ts/2);

sum_signal_in_clean = sum(signal_in,1);
noise_signals_in = sig*randn(1,length(codes_in));
sum_signal_in = sum_signal_in_clean + noise_signals_in;

% Initialize CCAF matrix
ccaf_clean = zeros(length(tau_corr),length(dopp_corr));
ccaf = ccaf_clean;

% Generate receiver replicas of code and carrier and CCAF
rep_code_prn = zeros(length(tau_corr),length(sum_signal_in));
for i = 1: length(tau_corr)
    rep_code_prn(i,:) = circshift(codes(sv_index_spoofed,:),i_tau_corr(i));
end
for j = 1 : length(dopp_corr)
    fprintf('%4d out of %4d\n',j,length(dopp_corr));
    rep_carr_prn = repmat(exp(-1i*2*pi*dopp_corr(j)*t),length(tau_corr),1);
    % OUTPUT
    ccaf_clean(:,j) = sum_signal_in_clean*transpose(rep_code_prn.*rep_carr_prn);   
    ccaf(:,j) = sum_signal_in*transpose(rep_code_prn.*rep_carr_prn);   
end

ccaf_clean = ccaf_clean/length(sum_signal_in);     
ccaf = ccaf/length(sum_signal_in);     

% Generate plots
% if (length(dopp_corr) > 1 && length(tau_corr) > 1)
% 
%     % Plot |CCAF|
%     figure(1)
%     surf(dopp_corr,chip_frac,abs(ccaf));shg
%     xlabel('Doppler (Hz)')
%     ylabel('tau (chips)')
%     title('|CCAF|')
% 
%     % Plot Re(CCAF)
%     figure(2)
%     surf(dopp_corr,chip_frac,abs(ccaf));shg
%     xlabel('Doppler (Hz)')
%     ylabel('tau (chips)')
%     title('Re(CCAF)')
% 
%     % Plot Im(CCAF)
%     figure(3)
%     surf(dopp_corr,chip_frac,imag(ccaf));shg
%     xlabel('Doppler (Hz)')
%     ylabel('tau (chips)')
%     title('Im(CCAF)')
% 
% end
