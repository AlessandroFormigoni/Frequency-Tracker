dt = 0.0001;
t = 0:dt:0.2;
fs = 1/dt;
fc = 200;
y = sin(2*pi*30*t)+2*sin(2*pi*60*t);
fDev = 50;
yMod = fmmod(y,fc,fs,fDev);

samples = length(y);
x = zeros(3, samples);
P = zeros(3, 3, samples);
z = @(x) x(1)*sin(x(3));
H = @(x) [sin(x(3)), 0, x(1)*cos(x(3))];
F = [1 0 0; 0 1 0; 0 dt 1];
yMeas = zeros(1, samples);
est = zeros(1, samples);
freq = zeros(1, samples);

%~

%Varianze dei rumori
q = 1e0;
r = 1e-2;
a = 1e-1;
lambda = 7;
I = lambda.*[dt*a 0 0; 0 dt*q (dt)^2*q/2;0 (dt)*q^2/2 (dt)^3*q/3];


%Stato iniziale
x(:,1) = [1; fc; 0];

%Matrice di covarianza
P(:,:,1) = eye(3);



for k = 1:samples-1
    %x(t|t-1)
    xhat_t = x(:,k);
   
    %Misura
    yMeas(k) = yMod(k) + wgn(1,1,10*log10(r));

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

end

figure(1)
plot(t,yMod,t,x(1,:))
figure(2)
ze = fmdemod(x(1,:),fc,fs,fDev);
plot(t,ze); hold on
zn = fmdemod(yMeas,fc,fs,fDev);
plot(t,zn)
z = fmdemod(yMod,fc,fs,fDev);
plot(t,z)


function sigmodel = f(x)
sigmodel = [2*cos(x(3))*x(1)-x(2); x(1); x(3)];
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