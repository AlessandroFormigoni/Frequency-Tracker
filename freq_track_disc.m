function [x_up, xhat_t, P_up, Kt, s, e_t] = freq_track_disc(yMeas,x,P,r,q,w)
H = [1 0 0];
I = diag([w w q]);
%~
    %errore
    e_t = yMeas - H*(x);

    %Guadagno del filtro
    s = (H*P*H' + r);
    Kt = P*H' / s;
    
    %Correzione
    xhat_t = x + Kt*e_t;

    %Measurement update
    F = linearize(xhat_t);
    L = chol(F*(P-Kt*H*P)*F'+ I);
    P_up = L'*L;
    x_up = f(xhat_t);
end

function sigmodel = f(x)
sigmodel = [cos(x(3))*x(1)-sin(x(3))*x(2); sin(x(3))*x(1)+cos(x(3))*x(2); x(3)];
end

function Ft = linearize(xhat)
Ft = zeros(3,3);
Ft(1,:) = [cos(xhat(3)) -sin(xhat(3)) (-sin(xhat(3))*xhat(1)-cos(xhat(3))*xhat(2))];
Ft(2,:) = [sin(xhat(3)) cos(xhat(3)) (cos(xhat(3))*xhat(1)-sin(xhat(3))*xhat(2))];
Ft(3,:) = [0 0 1];
end