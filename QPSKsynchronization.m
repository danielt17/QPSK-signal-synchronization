%% An educational code about syncronization based on the guide in PySDR, written completely in matlab
% https://pysdr.org/content/sync.html
% Adapted by Daniel Teitelman

%% Setup

close all;
clear;
clc;

%% Signal generation 

% Generating QPSK waveform

num_symbols = 200;
sps = 8;
symbols = randsrc(1,num_symbols,[1+1j 1-1j -1+1j -1-1j])/sqrt(2);
pulse_train = [];
for i  = 1:num_symbols
    pulse = zeros(1,sps);
    pulse(1) = symbols(i);
    pulse_train = [pulse_train,pulse];
end

% Passing the pulse train through a raised cosine filter - this is done to 
% minimize intersymbol interference (ISI). 
% We first create the filter than convolve it with the pulse train, thus we
% pulse shaped our waveform.

% Create pulse shaping filter

num_taps = 101;
beta = 0.35;
Ts = sps; % Assume sample rate is 1 Hz, so sample period is 1, so *symbol* period is 8
t = -51:1:51;
h = sinc(t./Ts) .* cos(pi.*beta.*t./Ts) ./ (1 - (2.*beta.*t./Ts).^2);

% Filter the signal

samples = conv(pulse_train,h,'same');

% Plot waveform generation pipeline

figure;
subplot(3, 2, 1)
plot(real(symbols),imag(symbols),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('Transmitted QPSK symbols','FontSize',12)
subplot(3, 2, 2)
plot(real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('Pulsed shaped signal','FontSize',12)
subplot(3, 2, 3)
plot(0:1:num_symbols-1,real(symbols),0:1:num_symbols-1,imag(symbols),'LineWidth',2)
grid on;
xlabel('Time [sec]','FontSize',12)
ylabel('Amplitude [U.O.S]','FontSize',12)
title('Transmitted QPSK symbols. I\\Q vs time plot','FontSize',12)
legend('In-Phase','Quadrature')
subplot(3, 2, 4)
plot(get_t(samples),real(samples),get_t(samples),imag(samples),'LineWidth',2)
grid on;
xlabel('Time [sec]','FontSize',12)
ylabel('Amplitude [U.O.S]','FontSize',12)
title('Pulsed shaped signal. I\\Q vs time plot','FontSize',12)
legend('In-Phase','Quadrature')
subplot(3, 2, [5 6])
stem(t,h,'LineWidth',2)
grid on;
xlabel('Taps [#]','FontSize',12)
ylabel('Amplitude','FontSize',12)
title('Raised cosine filter taps','FontSize',12)

%% Simulating the channel

% We add a delay, an operation which is done by using a fractional delay
% filter. For further reading look no further than https://www.intechopen.com/chapters/18566

delay = 0.4; % fractional delay, in samples
N = 21; % number of tps
n = floor(-N/2):floor(N/2)-1; % ...-3,-2,-1,0,1,2,3...
h = sinc(n-delay); % calc filter taps
h = h .* hamming(N)'; % # window the filter to make sure it decays to 0 on both sides
h = h./ sum(h); % normalize to get unity gain, we don't want to change the amplitude/power
samples_old = samples;
samples = conv(samples,h,'same'); % apply filter

% Plot delay channel result

figure;
subplot(4, 2, [1 2 3 4])
plot(real(samples_old),imag(samples_old),'o',real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('I\\Q plot - Comparison between before and after applying the fractional delay filter','FontSize',12)
legend('Before fractional delay filter','After fractional delay filter')
subplot(4, 2, 5)
plot(get_t(samples_old(1:200)),real(samples_old(1:200)),get_t(samples(1:200)),real(samples(1:200)),'LineWidth',2)
grid on;
xlabel('Time [sec]','FontSize',12)
ylabel('Amplitude [U.O.S]','FontSize',12)
title('In-Phase vs time','FontSize',12)
legend('Before fractional delay filter','After fractional delay filter')
subplot(4, 2, 6)
plot(get_t(samples_old(1:200)),imag(samples_old(1:200)),get_t(samples(1:200)),imag(samples(1:200)),'LineWidth',2)
grid on;
xlabel('Time [sec]','FontSize',12)
ylabel('Amplitude [U.O.S]','FontSize',12)
title('Quadrature vs time','FontSize',12)
legend('Before fractional delay filter','After fractional delay filter')
subplot(4, 2, [7 8])
stem(n,h,'LineWidth',2)
grid on;
xlabel('Taps [#]','FontSize',12)
ylabel('Amplitude','FontSize',12)
title('Fractional delay filter taps','FontSize',12)

% Adding a frequency offset

fs = 1e6; % assume our sample rate is 1 MHz
fo = 13000; % simulte freq offset
Ts = 1/fs; % calc sample period
t = 0:Ts:Ts*(size(samples,2)-1);
samples_old = samples;
samples = samples .* exp(1j.*2.*pi.*fo.*t);

% Plot frequency offset channel result

figure;
subplot(3, 2, [1 2 3 4])
plot(real(samples_old),imag(samples_old),'o',real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('I\\Q plot - Comparison between before and after applying a frequency offset channel','FontSize',12)
legend('Before frequency offset channel','After frequency offset channel')
subplot(3, 2, 5)
plot(get_t(samples_old(1:600)),real(samples_old(1:600)),get_t(samples(1:600)),real(samples(1:600)),'LineWidth',2)
grid on;
xlabel('Time [sec]','FontSize',12)
ylabel('Amplitude [U.O.S]','FontSize',12)
title('In-Phase vs time','FontSize',12)
legend('Before frequency offset channel','After frequency offset channel')
subplot(3, 2, 6)
plot(get_t(samples_old(1:600)),imag(samples_old(1:600)),get_t(samples(1:600)),imag(samples(1:600)),'LineWidth',2)
grid on;
xlabel('Time [sec]','FontSize',12)
ylabel('Amplitude [U.O.S]','FontSize',12)
title('Quadrature vs time','FontSize',12)
legend('Before frequency offset channel','After frequency offset channel')

% AWGN channel

samples_old = samples;
samples = awgn(samples,15,'measured');

figure;
plot(real(samples_old),imag(samples_old),'o',real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('I\\Q plot - Comparison between before and after applying an AWGN channel','FontSize',12)
legend('Before AWGN channel','After AWGN channel')

%% Time synchronization - Based on Mueller and Muller Timing Synchronization Algorithm
% further reading https://wirelesspi.com/mueller-and-muller-timing-synchronization-algorithm/

% Time synchronization is done to deal with the fractional time delay
% filter, we start with interpolating the signal, and afterwards apply the
% timing synchronization algorithm

samples_interpolated = resample(samples, 16, 1);

mu = 0; % initial estimate of phase of the sample
out = zeros(1,length(samples)+10);
out_rail = zeros(1,length(samples)+10); % stores values, each iteration we need the previous 2 values plus current value
i_in = 1; % input samples index
i_out = 3; % output index (let first two outputs be 0)
while (i_out < length(samples)) && (i_in + 16 < length(samples))
    out(i_out) = samples_interpolated(i_in*16 + round(mu*16)); % grab what we think is the "best" sample
    out_rail(i_out) = round(real(out(i_out)) > 0) + 1j*round(imag(out(i_out)) > 0);
    x = (out_rail(i_out) - out_rail(i_out-2)) * conj(out(i_out-1));
    y = (out(i_out) - out(i_out-2)) * conj(out_rail(i_out-1));
    mm_val = real(y-x);
    mu = mu + sps + 0.3*mm_val;
    i_in = i_in + floor(mu); % round down to nearest int since we are using it as an index
    mu = mu - floor(mu); % remove the intger part of mu
    i_out = i_out + 1;
end
out = out(3:i_out+1);
samples_old = samples;
samples = out;

% Plot timing synchronization algorithm result

figure;
plot(real(samples_old),imag(samples_old),'o',real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('I\\Q plot - Comparison between before and after applying Mueller and Muller timing synchronization algorithm','FontSize',12)
legend('Before timing synchronization','After timing synchronization')

%% Coarse frequency synchronization

% We start with coarse frequency synchronization, we should remember  our
% signal is built from QPSK symbols, therefore one could look at the
% following simple trick to get a rough estimate of the frequency offset 
% the signal experienced.
% r(t) = s(t)*e^(2*pi*fo*t) + n(t) if we do two sequences of complex
% multiplication we get r^4(t) = e^(2*pi*4*fo*t) as s(t)^4 = 1.

freqs = -fs/2:fs/2;
fft_samples = abs(fftshift(fft(samples,length(freqs))));
fft_samples_squared = abs(fftshift(fft(samples .^ 2,length(freqs))));
fft_samples_quad  = abs(fftshift(fft(samples .^ 4,length(freqs))));
[max_freq_fft_value,max_freq_ind] = max(fft_samples_quad);
max_freq = freqs(max_freq_ind);
frequency_offset = max_freq/4;

figure;
subplot(4, 1, 1)
plot(freqs*10^(-6),fft_samples,'LineWidth',2)
grid on;
xlabel('Frequency [Mhz]','FontSize',12)
ylabel('Amplitude','FontSize',12)
title('PSD of received signal','FontSize',12)
subplot(4, 1, 2)
plot(freqs*10^(-6),fft_samples_squared,'LineWidth',2)
grid on;
xlabel('Frequency [Mhz]','FontSize',12)
ylabel('Amplitude','FontSize',12)
title('PSD of received signal squared','FontSize',12)
subplot(4, 1, 3)
plot(freqs*10^(-6),fft_samples_quad,'LineWidth',2); hold on;
stem(max_freq*10^(-6),max_freq_fft_value,'LineWidth',1.5)
legend('PSD','Frequency shift estimator','Location','northwest')
grid on;
xlabel('Frequency [Mhz]','FontSize',12)
ylabel('Amplitude','FontSize',12)
title('PSD of received signal quadropled','FontSize',12)
subplot(4, 1, 4)
plot(freqs(max_freq_ind-5e4:max_freq_ind+5e4)*10^(-6),fft_samples_quad(max_freq_ind-5e4:max_freq_ind+5e4),'LineWidth',2); hold on;
stem(max_freq*10^(-6),max_freq_fft_value,'LineWidth',1.5)
legend('PSD','Frequency shift estimator','Location','northwest')
grid on;
xlabel('Frequency [Mhz]','FontSize',12)
ylabel('Amplitude','FontSize',12)
title('PSD of received signal quadropled - Zoom in','FontSize',12)

Ts = 1/fs;
t = 0:Ts:Ts*(size(samples,2)-1);
samples_old = samples;
samples = samples.*exp(-1j*2*pi*frequency_offset*t);

figure;
plot(real(samples_old),imag(samples_old),'o',real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('I\\Q plot - Comparison between before and after applying coarse frequency synchronization','FontSize',12)
legend('Before coarse frequency synchronization','After coarse frequency synchronization')


%% Fine Frequency Synchronization

% Fine frequency synchronization is done in order to minimize the small
% error that we didn't deal with using the previous methods. We use a
% feedback loop called a Costas Loop which is a type of PLL, for further
% reading one should look at the following source: https://en.wikipedia.org/wiki/Costas_loop

N = length(samples);
phase = 0;
freq = 0;
% These next two params is what to adjust, to make the feedback loop faster or slower (which impacts stability)
alpha = 0.132;
beta = 0.00932;
out = zeros(1,N);
freq_log = [];
for i = 1:N
    out(i) = samples(i) * exp(-1j*phase); % adjust the input sample by the inverse of the estimated phase offset
    error = order_4_costas_loop_error_equation(out(i)); % This is the error formula for 2nd order Costas Loop (e.g. for BPSK)
    % Advance the loop (recalc phase and freq offset)
    freq = freq + beta*error;
    freq_log = [freq_log; (freq*fs/(2*pi))]; % convert from angular velocity to Hz for logging
    phase = phase + freq + (alpha * error);
    while phase >= 2*pi
        phase = phase - 2*pi;
    end
    while phase < 0
        phase = phase + 2*pi;
    end
end

samples_old = samples;
samples = out;

figure;
subplot(4,2,[1 2 3 4 5 6])
plot(real(samples_old),imag(samples_old),'o',real(samples),imag(samples),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title('I\\Q plot - Comparison between before and after applying fine frequency synchronization','FontSize',12)
legend('Before fine frequency synchronization','After fine frequency synchronization')
subplot(4,2,[7 8])
plot(freq_log,'-','LineWidth',2)
grid on;
xlabel('Sample [#]','FontSize',12)
ylabel('Frequency offset [Hz]','FontSize',12)
title('Frequency offset with respect to sample number [#]','FontSize',12)

%% Evaluating the result

% We estimate the EVM of the synchronization procedure, as defined here 
% https://en.wikipedia.org/wiki/Error_vector_magnitude and compare the
% result with our starting symbols

EVM_value = EVM(samples,symbols);

figure;
plot(real(samples),imag(samples),'o',real(symbols),imag(symbols),'o','LineWidth',2)
grid on;
xlabel('In-Phase','FontSize',12)
ylabel('Quadrature','FontSize',12)
title(strcat('I\\Q plot - Comparison between the transmitted signals and the final synchronized signal at the receiver. EVM = ',num2str(EVM_value),' [dB]'),'FontSize',12)
legend('Received synchronized signal','Transmitted signal')


%% Functions

function t = get_t(y)
% create t linspace for given size
    t = 0:1:(size(y,2)-1);
end

function error = order_4_costas_loop_error_equation(sample)
    if real(sample) > 0
        a = 1;
    else
        a = -1;
    end
    if imag(sample) > 0
        b = 1;
    else
        b = -1;
    end
    error = a * imag(sample)  - b * real(sample);
end

function EVM_value = EVM(s_hat,s)
    EVM_value = 10 * log10(mean(abs(s_hat-s).^2)); 
end