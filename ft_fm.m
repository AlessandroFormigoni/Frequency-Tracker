clear
clc

dt = 0.0001;                 %Periodo di campionamento
t = (0:dt:2);                %Asse dei tempi
nsampl = length(t);         %Numero di campioni
x = zeros(3, nsampl);       %Vettore degli stati
P = zeros(3,3,nsampl);      %Matrice di covarianza

%Generazione segnale
fs = 1/dt;
fc = 200;
y = sin(2*pi*30*t)+2*sin(2*pi*60*t);
fDev = 50;
yMod = fmmod(y,fc,fs,fDev);
yMeas = zeros(1,nsampl);

%Condizioni iniziali
x(:,1) = [1 0 fDev];
P(:,:,1) = eye(3);

%Inizializzazione parametri
lambda = 5e-4;               %Parametro regolabile per la performance del filtro
r = 1e-2;                    %Varianza del rumore sulla misura
%v = 2*pinknoise(nsampl);    %Generazione rumore
v = wgn(1,nsampl,10*log10(r));
%Trasmissione segnale su canale rumoroso
yTrans = yMod + v;
z = fmdemod(yTrans,fc,fs,fDev);
%Esecuzione algoritmo
for k = 1:nsampl-1
yMeas(k) = z(k);
[x(:,k+1), P(:,:,k+1)] = freq_track(dt,yMeas(k),x(:,k),P(:,:,k),lambda);
end

%Grafici
figure(1)
%plot(t,yMod,t,yMeas);
%hold on
plot(t,x(1,:),'k','LineWidth',2);
hold on
plot(t,y);
plot(t,z);
legend('y stimata','y');
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
varf = zeros(1,nsampl);
for k=1:nsampl
    varf(k) = P(3,3,k);
end
plot(t,varf);

sprintf('Errore quadratico medio della stima del segnale: %f', mean((y-x(1,:)).^2)/(var(y)))