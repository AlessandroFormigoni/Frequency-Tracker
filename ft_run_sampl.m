clear
clc
% nmsefreq = zeros(1,1000);
% nmsesig = zeros(1,1000);
passed = 0;
rmsefreq = zeros(120,200);
% cj = 1;
% for j = 1:120
% ci = 1;
%   for i=0.01:0.01:2
    dt = 1/(100);               %Periodo di campionamento
    t = 0:dt:16;                %Asse dei tempi
    nsampl = length(t);         %Numero di campioni
    x = zeros(3, nsampl);       %Vettore degli stati
    P = zeros(3,3,nsampl);      %Matrice di covarianza
    K = zeros(3,nsampl-1);
    s = zeros(1,nsampl-1);
    e = zeros(1,nsampl-1);
    scale_factor = 3;
    %Frequenza reale (iniziale)
    omega = 2*pi*1;
    err = 0.4*rand() - 0.2;              %Errore della frequenza sulla stima iniziale in %
    %Generazione segnale
    y = (t>=0 & t<=8).*cos(omega*t) + (t>8 & t<=16).*cos(scale_factor*omega*t) + (t>16).*cos(omega*t);
    
    yMeas = zeros(1,nsampl); 
    
    %Condizioni iniziali
    x(:,1) = [1 0 omega-err*omega];
    P(:,:,1) = eye(3);
    
    %Inizializzazione parametri
    q = 1;               %Varianza del rumore sulla frequenza
    r = 1e-1;            %Varianza del rumore sulla misura
    w = 1e-4;            %Varianza del rumore sull'ampiezza
    %Generazione rumore
    %v = 4*pinknoise(nsampl);
    v = wgn(1,nsampl,10*log10(r));
    %Esecuzione algoritmo
    for k = 1:nsampl-1
    yMeas(k) = y(k) + v(k);
    [x(:,k+1),x(:,k), P(:,:,k+1), K(:,k), s(k), e(k)] = freq_track(dt,yMeas(k),x(:,k),P(:,:,k),r,q,w);
    end
%   realFreq = (t>=0 & t<=8).*((omega)*ones(1, nsampl)) + (t>8).*((scale_factor*omega)*ones(1, nsampl));
    %Grafici
     figure(1)
      omega = omega/(2*pi);
      x(3,:) = x(3,:)./(2*pi);
      realFreq = (t>=0 & t<=8).*((omega)*ones(1, nsampl)) + (t>8 & t<=16).*((scale_factor*omega)*ones(1, nsampl)) + (t>16).*((omega)*ones(1, nsampl));
   subplot(2,1,1)
    plot(t,y);
    hold on
    plot(t,yMeas);
    plot(t,x(1,:),'k','LineWidth',2);
    legend('y','y misurata','y stimata');
    grid on
    subplot(2,1,2)
    plot(t,realFreq,'--',t,x(3,:));
    legend('frequenza reale', 'frequenza stimata');
    grid on
    
    figure(2)
    subplot(2,1,1)
    a = abs(x(1,:)+1i*x(2,:));
    plot(t, a)
    xlabel('t')
    ylabel('Ampiezza');
    grid on
    subplot(2,1,2)
    arg = angle(x(1,:)+1i*x(2,:));
    plot(t,arg);
    grid on;
    xlabel('t')
    ylabel('Fase');
    
    figure(3)
    subplot(2,1,1)
    tr = zeros(1,nsampl);
    for k=1:nsampl
        tr(k) = trace(P(:,:,k));
    end
    plot(t,tr);
    xlabel('t');
    ylabel('tr(P_t)')
    
    subplot(2,1,2)
    [c, lags] = xcorr(e,'coeff');
    plot(dt.*lags,c);
    xlabel('t');
    ylabel('E[e(t)e(t+\tau)]')

    %test std e chi^2
     E = e.*s.^(-1).*e;
     mE = nsampl*mean(E);
     r1 = chi2inv(0.025, nsampl);
     r2 = chi2inv(1-0.025, nsampl);
    figure(4)
    plot(e);
    xlabel('campioni')
    hold on
    grid on
    plot(2.*sqrt(s),'r')
    plot(-2.*sqrt(s),'r')
    ninf = sum(e<-2.*sqrt(s));
    sinf = sum(e>2.*sqrt(s));
    sprintf('Il %f ha passato il test std', (1-(ninf+sinf)/nsampl)*100)
   legend('e(n)','\pm 2 \surd{S(n)}')
    if mE <= r2 && mE >= r1
        sprintf('Il filtro ha passato il test chiquadro')
        passed = passed + 1;
    else
        sprintf('Il filtro non ha passato il test chiquadro')
    end
    
%     nmsefreq(i) = mean((realFreq-x(3,:)).^2)/(mean(realFreq.^2));
%      nmsesig(i) = mean((y-x(1,:)).^2)/(mean(y.^2));
%     rmsefreq(cj,ci) = sqrt(mean((realFreq-x(3,:)).^2));
    sprintf("Lambda = %2f", r/q)
    sprintf('NMSE della stima della frequenza: %f', mean((realFreq-x(3,:)).^2)/(mean(realFreq.^2)))
    sprintf('NMSE della stima del segnale: %f', mean((y-x(1,:)).^2)/(mean(y.^2)))
    sprintf('La media dell innovazione: %f', mean(e))
%     ci = ci + 1;
%   end
%   cj = cj + 1;
% end

% a = mean(nmsefreq);
% b = std(nmsefreq);
% c = mean(nmsesig);
% d = std(nmsesig);
