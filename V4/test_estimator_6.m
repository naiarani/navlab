%%  Reproduce true signal

St = true_amplitude*(1 - abs(true_chip_delay - chip_frac))';  St(St<0) = 0;
St_all = repelem(St,2*ldopp);

Sfp1 = sinc((dopp_corr - true_dopp_offset)*Tcohs)';
Sfp2 = exp(1i*(pi*(dopp_corr - true_dopp_offset)*Tcohs - true_phase))';
Sfp = Sfp1.*Sfp2;
Sfp_big = repmat(Sfp,length(chip_frac),1);
re_Sfp_big= reshape(real(Sfp_big)',[],1);
im_Sfp_big = reshape(imag(Sfp_big)',[],1);
Sfp_all = zeros(N,1); 
Sfp_all(1:2:end,1) = re_Sfp_big;
Sfp_all(2:2:end,1) = im_Sfp_big;

z_true = St_all.*Sfp_all;
% figure; plot(z_true,'.-'); shg

%%  True signal H

% Tau
dt = 1e-3;
St_ = true_amplitude*(1 - abs(true_chip_delay - (chip_frac - dt)))';  St_(St_<0) = 0;
dStdt = (St - St_)/dt;
dStdt_all = repelem(dStdt,2*ldopp);

% Doppler
df = 1e-3;
Sfp1_ = sinc(((dopp_corr  - df) - true_dopp_offset)*Tcohs)';
dSfp1df = (Sfp1 - Sfp1_)/df;
dSfp2df = 1i*pi*Tcohs*Sfp2;
dSfpdf = dSfp1df.*Sfp2 + Sfp1.*dSfp2df;
dSfpdf_big = repmat(dSfpdf,length(chip_frac),1);
re_dSfpdf_big = reshape(real(dSfpdf_big)',[],1);
im_dSfpdf_big = reshape(imag(dSfpdf_big)',[],1);
dSfpdf_all = zeros(N,1); 
dSfpdf_all(1:2:end,1) = re_dSfpdf_big;
dSfpdf_all(2:2:end,1) = im_dSfpdf_big;

% Phase
dSfp1dp = 0;
dSfp2dp = -1i*Sfp2;
dSfpdp = dSfp1dp.*Sfp2 + Sfp1.*dSfp2dp;
dSfpdp_big = repmat(dSfpdp,length(chip_frac),1);
re_dSfpdp_big = reshape(real(dSfpdp_big)',[],1);
im_dSfpdp_big = reshape(imag(dSfpdp_big)',[],1);
dSfpdp_all = zeros(N,1); 
dSfpdp_all(1:2:end,1) = re_dSfpdp_big;
dSfpdp_all(2:2:end,1) = im_dSfpdp_big;

H1 = St_all.*Sfp_all/true_amplitude;
H2 = dStdt_all.*Sfp_all;
H3 = dSfpdf_all.*St_all;
H4 = dSfpdp_all.*St_all;

H = [H1 H2 H3 H4];

%%  Reproduce spoofed signal

St = spoofed_amplitude*(1 - abs(spoofed_chip_delay - chip_frac))';  St(St<0) = 0;
St_all = repelem(St,2*ldopp);

Sfp1 = sinc((dopp_corr - spoofed_dopp_offset)*Tcohs)';
Sfp2 = exp(1i*(pi*(dopp_corr - spoofed_dopp_offset)*Tcohs - spoofed_phase))';
Sfp = Sfp1.*Sfp2;
Sfp_big = repmat(Sfp,length(chip_frac),1);
re_Sfp_big= reshape(real(Sfp_big)',[],1);
im_Sfp_big = reshape(imag(Sfp_big)',[],1);
Sfp_all = zeros(N,1); 
Sfp_all(1:2:end,1) = re_Sfp_big;
Sfp_all(2:2:end,1) = im_Sfp_big;

z_spoofed = St_all.*Sfp_all;
% figure; plot(z_spoofed,'.-'); shg

%%  Spoofed signal H

% Tau
dt = 1e-3;
St_ = spoofed_amplitude*(1 - abs(spoofed_chip_delay - (chip_frac - dt)))';  St_(St_<0) = 0;
dStdt = (St - St_)/dt;
dStdt_all = repelem(dStdt,2*ldopp);

% Doppler
df = 0.001;
dSfp1df = (Sfp1 - sinc(((dopp_corr - df) - spoofed_dopp_offset)*Tcohs)')/df;
dSfp2df = 1i*pi*Tcohs*Sfp2;
dSfpdf = dSfp1df.*Sfp2 + Sfp1.*dSfp2df;
dSfpdf_big = repmat(dSfpdf,length(chip_frac),1);
re_dSfpdf_big = reshape(real(dSfpdf_big)',[],1);
im_dSfpdf_big = reshape(imag(dSfpdf_big)',[],1);
dSfpdf_all = zeros(N,1); 
dSfpdf_all(1:2:end,1) = re_dSfpdf_big;
dSfpdf_all(2:2:end,1) = im_dSfpdf_big;

% Phase
dSfp1dp = 0;
dSfp2dp = -1i*Sfp2;
dSfpdp = dSfp1dp.*Sfp2 + Sfp1.*dSfp2dp;
dSfpdp_big = repmat(dSfpdp,length(chip_frac),1);
re_dSfpdp_big = reshape(real(dSfpdp_big)',[],1);
im_dSfpdp_big = reshape(imag(dSfpdp_big)',[],1);
dSfpdp_all = zeros(N,1); 
dSfpdp_all(1:2:end,1) = re_dSfpdp_big;
dSfpdp_all(2:2:end,1) = im_dSfpdp_big;

H1 = St_all.*Sfp_all/spoofed_amplitude;
H2 = dStdt_all.*Sfp_all;
H3 = dSfpdf_all.*St_all;
H4 = dSfpdp_all.*St_all;

H = [H H1 H2 H3 H4];                     %  Comment if no spoofing

z1 = z - (z_true + z_spoofed) ;       %  Comment if no spoofing
% z1 = z - z_true;                           %  Uncomment if no spoofing
% figure(7); stem(z1,'.'); shg; grid       % Plot un-whitened  measurements       

%% Estimator
  
% Without decorrelation
Hpi = pinv(H);
P1 = Hpi*COV*Hpi';
eP1 = eig(P1);
xhat1 = (H\z1)'
sigdiag1 = sqrt(diag(P1))'
ratio1 = abs(xhat1)./sqrt(diag(P1))'
r1 = norm(z1 - H*xhat1');
x1_norm = norm(xhat1);
v1_norm = norm(sigdiag1);

% With decorrelatiion
z2 = sqrtmCOV\z1;
% figure(8); stem(z2,'.'); shg; grid       % Plot whitened measurements            
Hw = sqrtmCOV\H;
P2 = inv(Hw'*Hw);
eP2 = eig(P2);
xhat2 = (Hw\z2)'
sigdiag2 = sqrt(diag(P2))'
ratio2 = abs(xhat2)./sqrt(diag(P2))'
r2 = norm(z2 - Hw*xhat2');
x2_norm = norm(xhat2);
v2_norm = norm(sigdiag2);
% hw_norm = norm(eye(length(z2)) - Hw*inv((Hw')*Hw)*Hw');
