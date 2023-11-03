clear
clc

filename = "fick_yer_mini_dbass.wav";
[y, Fs] = audioread(filename);
dt = 2*pi/Fs;
nsampl = length(y(:,1));
x = zeros(3,nsampl);
P = zeros(3,3,nsampl);
yMeasl = zeros(1,nsampl);
yMeasr = zeros(1,nsampl);
omega = 110;

xl(:,1) = [0 0 omega];
Pl(:,:,1) = eye(3);
xr(:,1) = [0 0 omega];
Pr(:,:,1) = eye(3);

r = 1e-2; lambda = 0.003;

for k = 1:nsampl-1
yMeasl(k) = y(k,1) + r*sin(100/201*k);%wgn(1,1,10*log10(r));
yMeasr(k) = y(k,2) + wgn(1,1,10*log10(r));
[xl(:,k+1), Pl(:,:,k+1)] = freq_track(dt,yMeasl(k),xl(:,k),Pl(:,:,k),lambda);
[xr(:,k+1), Pr(:,:,k+1)] = freq_track(dt,yMeasr(k),xr(:,k),Pr(:,:,k),lambda);
end

figure(1)
plot((xl(3,:)))
grid on
figure(2)
hold on
plot(y);
xs = [xl(1,:); xr(1,:)];
%plot(xs);
%grid on
%sound(yMeas,Fs);
%pause(10)
sound(xs, Fs);
