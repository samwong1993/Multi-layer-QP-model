%% Calculate A/B/C for F Layer
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [A B C] = ABC(F,R,Rb,Rm,Ym,beta)
A = 1-1/F^2+(Rb/F/Ym)^2;
B = - 2*Rm*Rb^2/(F^2*Ym^2);
C = (Rb*Rm/F/Ym)^2-R^2.*cos(beta).^2;
B = B*ones(1,length(C));
A = A*ones(1,length(C));
end