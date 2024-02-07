clear
clc

filename = "./music/1000.mp4";
[y, Fs] = audioread(filename);
dt = 1/(Fs);
nsampl = length(y(:,1));
x = zeros(3,nsampl);
P = zeros(3,3,nsampl);
yMeasl = zeros(1,nsampl);
yMeasr = zeros(1,nsampl);
s = zeros(1,nsampl-1);
e = zeros(1,nsampl-1);
omega = 67*(2*pi);
t = 0:dt:(dt*nsampl);
xl(:,1) = [0 0 omega];
Pl(:,:,1) = eye(3);
% xr(:,1) = [0 0 omega];
% Pr(:,:,1) = eye(3);
r = 0.011; lambda = 800; w = 7e-4;
%vl = 4*pinknoise(nsampl);
vr = wgn(1,nsampl,10*log10(r));
for k = 1:nsampl-1
yMeasl(k) = y(k,1) + 0*vr(k);
% yMeasr(k) = y(k,2) + vr(k);
[xl(:,k+1),xl(:,k), Pl(:,:,k+1),~,s(k),e(k)] = freq_track(dt,yMeasl(k),xl(:,k),Pl(:,:,k),r,lambda,w);
% [xr(:,k+1),xr(:,k), Pr(:,:,k+1)] = freq_track(dt,yMeasr(k),xr(:,k),Pr(:,:,k),r,lambda);
end

figure(1)
plot(t(1:length(t)-1),(xl(3,:)./(2*pi)))
xlabel('t')
ylabel('Hz')
legend('frequenza stimata','frequenza reale')
grid on
figure(2)
hold on
plot(y);
xs = [xl(1,:)];
figure(3)
pspectrum(y(:,1), Fs, 'spectrogram', 'FrequencyLimits', [20 1500], 'OverlapPercent',0,'Leakage',0.20,'MinThreshold',-45);
%plot(xs);
%grid on
%sound(yMeas,Fs);
%pause(10)
sound(xs, Fs);
