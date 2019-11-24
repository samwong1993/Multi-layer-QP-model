%% Calculate A/B/C£¨E + Joining layer£©
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [A B C] = QP_ABC(R,Rm,Rb,Ym,F,beta,joining)
M = length(beta);
if joining == true
    A = 1 - 1/F^2 - (Rb/F/Ym)^2;
    B = 2*Rm*Rb^2/F^2/Ym^2;
    C = -(Rb*Rm/F/Ym)^2 - R^2*cos(beta).^2;
    A = A*ones(1,M);
    B = B*ones(1,M);
else
	A = 1 - 1/F^2 + (Rb/F/Ym)^2;
    B = - 2*Rm*Rb^2/F^2/Ym^2;
    C = (Rb*Rm/F/Ym)^2 - R^2*cos(beta).^2;
    A = A*ones(1,M);
    B = B*ones(1,M);
end
end