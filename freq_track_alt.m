function x = freq_track_alt()
dt = 0.001;
t = 0:dt:12;
omega = pi/2;
y = (t>=0 & t<=5).*sin(omega*t) + (t>5).*sin(4*omega*t+pi/2);
%y = sin(omega*t);
samples = length(y);
x = zeros(3, samples);
P = zeros(3, 3, samples);
z = @(x) x(1)*sin(x(3));
H = @(x) [sin(x(3)), 0, x(1)*cos(x(3))];
F = [1 0 0; 0 1 0; 0 dt 1];
yMeas = zeros(1, samples);
est = zeros(1, samples);
freq = zeros(1, samples);
realFreq = (t>=0 & t<=5).*((omega)*ones(1, samples)) + (t>5).*((4*omega)*ones(1, samples));
%realFreq = ((omega)*ones(1, samples));

%~

%Varianze dei rumori
q = 1e-2;
r = 0.1;
a = 1e-3;
lambda = 2;
I = lambda.*[dt*a 0 0; 0 dt*q (dt)^2*q/2;0 (dt)*q^2/2 (dt)^3*q/3];


%Stato iniziale
x(:,1) = [1; omega-pi/4; 0];

%Matrice di covarianza
P(:,:,1) = eye(3);



for k = 1:samples
    %x(t|t-1)
    xhat_t = x(:,k);
   
    %Misura
    yMeas(k) = y(k) + wgn(1,1,10*log10(r));

    %errore
    e_t = yMeas(k) - z(xhat_t);

    %Guadagno del filtro
    Kt = P(:,:,k)*H(xhat_t)'.*inv(H(xhat_t)*P(:,:,k)*H(xhat_t)' + r);

    %Correzione
    x(:,k) = xhat_t + Kt.*e_t;

    %Measurement update
    P(:,:,k) = ((eye(3) - Kt*H(xhat_t))*P(:,:,k));
    
    %Predizione
    P(:,:,k+1)= (F*P(:,:,k)*F'+I);
    x(:,k+1) = F*(x(:,k));

    est(k) = x(1,k)*sin(x(3,k));
    freq(k) = x(2,k);

end
%finfreq = movmean(x(2,:),1/dt);
subplot(2,1,1)
plot(t,y,t,yMeas,t,est);
subplot(2,1,2)
plot(t,freq,t,realFreq);
end

function sigmodel = f(x)
sigmodel = [2*cos(x(3))*x(1)-x(2); x(1); x(3)+wgn(1,1, 10*log10(0.005))];
end

function Ft = linearize(xhat)
Ft = zeros(3,3);
Ft(1,:) = [2*cos(xhat(3)) -1 -2*sin(xhat(3))*xhat(1)];
Ft(2,:) = [1 0 0];
Ft(3,:) = [0 0 1];
end

function Kt = kalman_gain(P, H, r)
Kt = P*H'.*(H*P*H' + r).^(-1);
end