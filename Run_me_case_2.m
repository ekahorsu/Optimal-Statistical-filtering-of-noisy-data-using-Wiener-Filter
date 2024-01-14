clear;
close all;



%%% providing dn %%% 
filename = 'one_new.m4a';
[voice, Fs] = audioread(filename, [20000, 60000]);
dn = voice;
N = length(voice);
n = 1:N;


%%% providing vn %%%
vn=0.004*randn(1,N);   %% Here the factor 0.004 is chosen such that the noise amplitude is essentially a fraction of the aplitude of dn
v1(1)=vn(1);
v2(1)=vn(1);
for i=2:N
    v1(i)=0.8*v1(i-1)+vn(i);    %primary sensor noise
    v2(i)=-0.6*v2(i-1)+vn(i);   %secondary sensor noise
    %both noises are correlated
end
v1 = v1';
v2 = v2';

xn=dn+v1;       % Recieved Signal
%%%%% checking whether xn is stationary or not %%%%%

%%%% this particular part of the code is not ours, we have used a code that we got from Matlab Central %%% 

% perform WSS test and visualize the results
% Note: the signal is estimated as wide-sense stationary when it is
% simultaneously stationary about its mean, variance and autocovariance.
commandwindow
gamma = 0.9;
if isstationary(xn, gamma)
    disp('The signal is wide-sense stationary!')
else
    disp('The signal is not wide-sense stationary!') 
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting the input signal %%
figure()
plot(1:N,xn(1:N));
hold on;
plot(1:N,dn(1:N),'r');
xlabel('Time'); ylabel('Amplitude'); grid on; axis tight;
legend('Signal with Primary Sensor Noise','Desired Signal');

%% plotting the secondary noise%%
figure();
plot(1:N,v2(1:N));
title('Plot of Secondary Sensor Noise');

%% plotting the autocorrelation of the secondary noise %%
rv2=autocorr(v2(1:N));
stem(rv2,'xr');
xlabel('Index');ylabel('Amplhitude');grid on; axis tight;
title('Autocorrelation of Secondary Sensor Noise Signal');

%% plotting the cross-correlation of the input signal and the secondary noise %%
rxv2=xcorr(xn(1:N),v2(1:N));
figure();
stem(rxv2);
xlabel('Index');
title('Cross-correlation between Received Signal and Secondary Noise');grid on; axis tight;

% Order 8 Wiener Filter

order = 8;
Rv2 = zeros(order,order);
for i = 1: order
    count = i;
    for j = 1 : order
        Rv2(i,j) = rv2(count);
        if (j<i)
            count = count -1;
        else
            count = count +1;
        end
        
    end
end

wiener2 = Rv2\rxv2(1:order);

v1_est = conv(wiener2,v2);
d_est = xn - v1_est(1:N);

figure();
plot(d_est(1:N));
hold on
plot(dn(1:N),'r');
hold off;
title('Estimated Signal and Desired Signal for Order = 8');
xlabel('Time'); ylabel('Amplitude'); grid on; axis tight;
legend('Estimated Signal', 'Original Signal');
MSE_d_8 = mse(d_est,dn);
MSE_v1_8 = mse(v1_est,v1);
display(MSE_d_8);
display(MSE_v1_8);

% Order 64 Wiener Filter

order = 16;
Rv2 = zeros(order,order);
for i = 1: order
    count = i;
    for j = 1 : order
        Rv2(i,j) = rv2(count);
        if (j<i)
            count = count -1;
        else
            count = count +1;
        end
        
    end
end

wiener2 = Rv2\rxv2(1:order);

v1_est = conv(wiener2,v2);

d_est = xn - v1_est(1:N);

figure();
plot(d_est(1:N));
hold on
plot(dn(1:N),'r');
hold off;
title('Estimated Signal and Desired Signal for order 64');
xlabel('Time'); ylabel('Amplitude'); grid on; axis tight;
legend('Estimated Signal', 'Original Signal');
MSE_d_64 = mse(d_est,dn);
MSE_v1_64 = mse(v1_est,v1);
display(MSE_d_64);
display(MSE_v1_64);

%% comparing the frequency response of the estimated signal with the desired signal %%
f1 = fft(d_est);
f2 = fft(dn);
figure()
plot(f1);
hold on 
plot(f2, 'r');
hold off
title('Frequency response of Estimated Signal and Desired Signal for order 64');
xlabel('\omega'); ylabel('fourier_cooefficients'); grid on; axis tight;
legend('Estimated Signal', 'Original Signal');


