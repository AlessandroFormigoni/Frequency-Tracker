function [x_up, P_up] = freq_track(dt,yMeas,x,P,lambda)
H = [1 0 0];
I = diag([0 0 dt]);
%~
    %errore
    e_t = yMeas - H*x;

    %Guadagno del filtro
    Kt = P*H' / (H*P*H' + lambda);

    %Correzione
    xhat_t = x + Kt.*e_t;

    %Measurement update
    F = linearize(dt.*xhat_t);
    P_up = F*(P-Kt*H*P)*F'+ I;

    x_up = f([xhat_t(1) xhat_t(2) xhat_t(3)*dt]);
    x_up(3) = xhat_t(3);

end

function sigmodel = f(x)
sigmodel = [cos(x(3))*x(1)-sin(x(3))*x(2); sin(x(3))*x(1)+cos(x(3))*x(2); x(3)];
end

function Ft = linearize(xhat)
Ft = zeros(3,3);
Ft(1,:) = [cos(xhat(3)) -sin(xhat(3)) -sin(xhat(3))*xhat(1)-cos(xhat(3))*xhat(2)];
Ft(2,:) = [sin(xhat(3)) cos(xhat(3)) cos(xhat(3))*xhat(1)- sin(xhat(3))*xhat(2)];
Ft(3,:) = [0 0 1];
end
