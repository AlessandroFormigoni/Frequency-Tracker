clear
clc

dt = 1;                 %Periodo di campionamento
t = 0:dt:1600;                %Asse dei tempi
nsampl = length(t);         %Numero di campioni
x = zeros(3, nsampl);       %Vettore degli stati
P = zeros(3,3,nsampl);      %Matrice di covarianza
s = zeros(1,nsampl-1);
e = zeros(1,nsampl-1);

%Frequenza reale (iniziale)
omega = pi/8;
err = 0.2 - 0.4*rand();              %Errore della frequenza sulla stima iniziale in %
%Generazione segnale
y = (t>=0 & t<=800).*(cos(omega*t)) + (t>800).*cos(3*omega*t);
yMeas = zeros(1,nsampl);

%Condizioni iniziali
x(:,1) = [1 0 omega-err*omega];
P(:,:,1) = eye(3);
K = zeros(3,nsampl);

%Inizializzazione parametri
r = 5e-2;                 %Varianza del rumore sulla misura
q = 3e-3;               %Parametro regolabile per la performance del filtro
w = 9e-4;
%v = 2*pinknoise(nsampl);    %Generazione rumore
v = wgn(1,nsampl,10*log10(r));
%Esecuzione algoritmo
for k = 1:nsampl-1
yMeas(k) = y(k) + v(k);
[x(:,k+1),x(:,k), P(:,:,k+1),K(:,k), s(k), e(k)] = freq_track_disc(yMeas(k),x(:,k),P(:,:,k),r,q,w);
end

%Grafici
figure(1)
subplot(2,1,1)
realFreq = (t>=0 & t<=800).*((omega)*ones(1, nsampl)) + (t>800).*((3*omega)*ones(1, nsampl));
plot(t,y);
hold on
plot(t,yMeas);
plot(t,x(1,:),'k','LineWidth',2);
legend('y','y misurata','y stimata');
grid on
subplot(2,1,2)
plot(t,realFreq,'--',t,x(3,:));
legend('frequenza reale', 'frequenza stimata');
xlabel('campioni')
ylabel('rad/campione')
grid on

figure(2)
subplot(2,1,1)
a = abs(x(1,:)+1i*x(2,:));
plot(t, a)
ylabel('Ampiezza');
grid on
subplot(2,1,2)
arg = angle(x(1,:)+1i*x(2,:));
plot(t,arg);
ylabel('Fase');

figure(3)
subplot(2,1,1)
tr = zeros(1,nsampl);
for k=1:nsampl
    tr(k) = trace(P(:,:,k));
end
plot(t,tr);
xlabel('t (samples)');
ylabel('tr(P_t)')
subplot(2,1,2)
[c, lags] = xcorr(e,'coeff');
plot(dt.*lags,dt.*c);
xlabel('t (samples)');
ylabel('correlazione normalizzata')

figure(4)
plot(e);
hold on
grid on
plot(2.*sqrt(s))
plot(-2.*sqrt(s))
samp = nsampl-1;
E = e(1:samp).*s(1:samp).^(-1).*e(1:samp);
mE = length(1:samp)*mean(E);
r1 = chi2inv(0.025, samp);
r2 = chi2inv(1-0.025, samp);
if mE <= r2 && mE >= r1
    sprintf('Il filtro ha passato il test chiquadro')
else
    sprintf('Il filtro non ha passato il test chiquadro')
end

sprintf("Lambda = %2f", r/q)
sprintf('NMSE della stima della frequenza: %f', mean((realFreq-x(3,:)).^2)/(mean(realFreq.^2)))
sprintf('NMSE della stima del segnale: %f', mean((y-x(1,:)).^2)/(mean(y.^2)))
sprintf('La media dell innovazione: %f', mean(e))