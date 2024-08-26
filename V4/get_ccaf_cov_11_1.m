% Create measurement vector: re_im_ccaf = [ I(tau1,fd1); Q(tau1,fd1); I(tau1,fd2);
% Q(tau1,fd2); ... ; I(tau2,fd1); Q(tau2,fd1); I(tau2,fd2); Q(tau2,fd2) ;... ] 

re_ccaf = reshape(real(ccaf)',[],1);
im_ccaf = reshape(imag(ccaf)',[],1);
N = 2*length(im_ccaf)

re_im_ccaf = zeros(N,1);
re_im_ccaf(1:2:end,1) = re_ccaf;
re_im_ccaf(2:2:end,1) = im_ccaf;
z = re_im_ccaf;
% figure(4); stem(z,'.-'); shg; grid   % Plot raw measurements
% yax4 = get(gca,'Ylim');

% Create covariance matrix

% Code phase (chip_frac) contribution
ltau = length(chip_frac);
ldopp = length(dopp_corr);
t1mt2 = meshgrid(chip_frac)' - meshgrid(chip_frac);
R = 1 - abs(t1mt2);
R(R<0) = 0;
R_all = repelem(R,2*ldopp,2*ldopp);

% Doppler (dopp_corr) contribution
f1mf2 = meshgrid(dopp_corr)' - meshgrid(dopp_corr);
f1pf2 = meshgrid(dopp_corr) + meshgrid(dopp_corr)';                                               
scm = sinc(2*f1mf2*Tcohs);                                     
scp = sinc(2*f1pf2*Tcohs);
ssii = scm + scp;
ssqq = scm - scp;
ssqi = -(sinc(f1mf2*Tcohs).*sin(pi*f1mf2*Tcohs) + sinc(f1pf2*Tcohs).*sin(pi*f1pf2*Tcohs));

ssiq = ssqi'; 
dssii = repelem(ssii,2,2);
S_ii = dssii.*repmat([1 0 ; 0 0],ldopp,ldopp);
dssqq = repelem(ssqq,2,2);
S_qq = dssqq.*repmat([0 0 ; 0 1],ldopp,ldopp);
dssiq = repelem(ssiq,2,2);
S_iq = dssiq.*repmat([0 1 ; 0 0],ldopp,ldopp);
dssqi = repelem(ssqi,2,2);
S_qi = dssqi.*repmat([0 0 ; 1 0],ldopp,ldopp); 
S_tot = S_ii + S_qq + S_iq + S_qi;
S_all = repmat(S_tot,ltau,ltau);

% Combine code phase and Doppler conributions
var = N0/4/Tcohs;
COV = R_all.*S_all*var;
eigCOV = eig(COV);
COV = COV + max(eigCOV)*N*eye(N)*eps;
sqrtmCOV = sqrtm(COV);   
N1 = rank(COV)
[Us,Ss,Vs] = svd(COV);
sdSs = sqrt(diag(Ss));
% figure(5); semilogy(sdSs,'o-'); grid




