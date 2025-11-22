clear;
clc;
close all;


N  = 4024;          % better to use even N
dt = 1;             
fs = 1/dt;          
t  = (0:N-1)'*dt;

k  = (0:N-1)';
f  = k*fs/N;

fL = 0.1;          
fH = 0.10;           % still unused, but ok for MB definition
idxMB  = (f>=fL & f<=fH);
idxLow = (f<fL & f>0);

% ---- custom flat PSD level ----
S_flat = 3e-4;      % <-- YOUR custom PSD value (units of whatever PSD you want)

S_noise = S_flat*ones(N,1);                 % flat part = S_flat
S_noise(idxLow) = S_flat*(fL./f(idxLow));   % 1/f below MB, scaled consistently

% Avoid DC issue
S_noise(1) = S_noise(find(f>0,1,'first'));

% ---- Generate colored noise ----
rng(1);
phi  = 2*pi*rand(N,1);
Spec = sqrt(S_noise).*exp(1j*phi);

% Hermitian symmetry (assuming N even)
Spec(N/2+2:end) = conj(flipud(Spec(2:N/2)));
Spec(1)       = real(Spec(1));
Spec(N/2+1)   = real(Spec(N/2+1));

noise = real(ifft(Spec))*sqrt(N);



%% Plot noise spectra before and after filtering
Nfft = N;
[PSD_in, f_psd]  = pwelch(noise,  [], [], Nfft, fs);

figure;
loglog(f_psd, PSD_in,  'DisplayName','Raw noise'); hold on;
xline([fL fH], ':',  'DisplayName','MB limits', 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('PSD');
legend('Location','best');
grid on;

figure;
plot(t, noise)