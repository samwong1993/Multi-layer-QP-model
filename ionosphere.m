%% Simulation for ionosphere(E + F +Joining layer)
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [P D] = ionosphere(R,Rb,Rm,Ym,A,B,C,beta,joining)
M = length(beta);
%Rt(index) = -0.99999999999*(B(index)+sqrt(B(index).^2 - 4*A(index).*C(index)))/2./A(index);
if joining == true
	Rt = ones(1,M)*Rm;
    D = R^2*cos(beta)./sqrt(-C).*(arcsinD(A,B,C,Rt) - arcsinD(A,B,C,Rb));
	Xb = real(A.*Rb.^2 + B.*Rb + C);
    Xm = real(A.*Rt.^2 + B.*Rt + C);
    P = (sqrt(Xm)-sqrt(Xb))./A + B./2./sqrt(-A.^3).*(arcsinP(A,B,C,Rt) - arcsinP(A,B,C,Rb));
else
    Rt = ones(1,M)*Rm;
    index = (B.^2 - 4*A.*C>=0);
    Rt(index) = -(B(index)+sqrt(B(index).^2 - 4*A(index).*C(index)))/2./A(index);
    if Rm<Rt
        Rt = Rm;
    end
    D = - R^2*cos(beta)./sqrt(C).*log(abs(lnD(A,B,C,Rt))./lnD(A,B,C,Rb));
    Xb = abs(A.*Rb.^2 + B.*Rb + C);
    Xm = abs(A.*Rt.^2 + B.*Rt + C);
    P = ( sqrt(Xm) - sqrt(Xb) )./A + B/2./sqrt(A.^3).*log(lnP(A,B,C,Xb,Rb)./lnP(A,B,C,Xm,Rt));
end
end



