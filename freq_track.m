function [x_up,xhat_t, P_up, Kt, s, e_t] = freq_track(dt,yMeas,x,P,r,q,w)
H = [1 0 0];
I = diag([w*dt w*dt q*dt]);
%~
    %errore
    e_t = yMeas - H*x;

    %Guadagno del filtro
    s = (H*P*H' + r);
    Kt = P*H' / s;
    
    %Correzione
    xhat_t = x + Kt.*e_t;

    %Measurement update
    F = linearize(xhat_t,dt);
    L = chol(F*(P-Kt*H*P)*F'+ I);
    P_up = L'*L;

    x_up = f(xhat_t,dt);
end

function sigmodel = f(x,dt)
sigmodel = [cos(dt*x(3))*x(1)-sin(dt*x(3))*x(2); sin(dt*x(3))*x(1)+cos(dt*x(3))*x(2); x(3)];
end

function Ft = linearize(xhat,dt)
Ft = zeros(3,3);
Ft(1,:) = [cos(dt*xhat(3)) -sin(dt*xhat(3)) dt*(-sin(dt*xhat(3))*xhat(1)-cos(dt*xhat(3))*xhat(2))];
Ft(2,:) = [sin(dt*xhat(3)) cos(dt*xhat(3)) dt*(cos(dt*xhat(3))*xhat(1)-sin(dt*xhat(3))*xhat(2))];
Ft(3,:) = [0 0 1];
end
