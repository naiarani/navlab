%% Reproduce true signal

Rt = (1 - abs(true_chip_delay - chip_frac))';  Rt(Rt<0) = 0;
Rt_all = true_amplitude*repelem(Rt,2*ldopp);

Ft = sinc((dopp_corr - true_dopp_offset)*Tcohs)';
Pt = exp(1i*(pi*(dopp_corr - true_dopp_offset)*Tcohs - true_phase))';
FPt = Ft.*Pt;
FPt_all = repmat(FPt,length(chip_frac),1);

re_FPt= reshape(real(FPt_all)',[],1);
im_FPt = reshape(imag(FPt_all)',[],1);

re_im_FPt = zeros(N,1);
re_im_FPt(1:2:end,1) = re_FPt;
re_im_FPt(2:2:end,1) = im_FPt;
re_im_FPRt = re_im_FPt.*Rt_all;

%% Reproduce spoofed signal

Rs = (1 - abs(spoofed_chip_delay - chip_frac))';  Rs(Rs<0) = 0;
Rs_all = spoofed_amplitude*repelem(Rs,2*ldopp);

Fs = sinc((dopp_corr - spoofed_dopp_offset)*Tcohs)';
Ps = exp(1i*(pi*(dopp_corr - spoofed_dopp_offset)*Tcohs - spoofed_phase))';
FPs = Fs.*Ps;
FPs_all = repmat(FPs,length(chip_frac),1);

re_FPs = reshape(real(FPs_all)',[],1);
im_FPs = reshape(imag(FPs_all)',[],1);

re_im_FPs = zeros(N,1);
re_im_FPs(1:2:end,1) = re_FPs;
re_im_FPs(2:2:end,1) = im_FPs;
re_im_FPRs = re_im_FPs.*Rs_all;

%% Reproduced combined signal

re_im_FPRts = re_im_FPRt + re_im_FPRs;
figure(14); plot(re_im_FPRts,'.-'); grid; shg
set(gca,'Ylim',yax4);

figure(24); plot(re_im_ccaf-re_im_FPRts,'.-'); shg; grid
set(gca,'Ylim',yax4);

zts = sqrtmCOV\re_im_FPRts;
figure(16); plot(zts,'.-'); grid; shg
% set(gca,'Ylim',yax6);
z_cov = sqrtmCOV\z;  
%% Error in reproduced combined signal relative to true input (no multipath)
dz=norm(z_cov - zts);
dz2 = norm(zts);
% figure(26) ;plot(dz,'.-'); grid; shg
% yax6 = get(gca,'Ylim');
% set(gca,'Ylim',yax6);

dz1= norm(z - re_im_FPRts)/norm(re_im_FPRts);
% figure(26) ;plot(dz,'.-'); grid; shg
% set(gca,'Ylim',yax6);

dx = norm(z_cov)/norm(zts);
dx1 = norm(z)/norm(re_im_FPRts);

% zt = sqrtmCOV\re_im_FPRt;
% figure(36); plot(zt,'.-'); grid; shg
% set(gca,'Ylim',yax6);
%   
% zs = sqrtmCOV\re_im_FPRs;
% figure(46); plot(zs,'.-'); grid; shg
% set(gca,'Ylim',yax6);
%   
% dzts=zt - zs;
% figure(56) ;plot(dzts,'.-'); grid; shg
% set(gca,'Ylim',yax6);