%% Calculate whether the penetrate the layer
%created by Huang Sen
%Email: huangsen1993@gmail.com
%joining  = 1 joining layer/ joining = 0 E/F layer
function [Upper] = penetrate(Rm,Rb,F,Ym,R,joining)
M = 1;
if joining == true
    if F<1
        con_beta = pi/2;
    else
        A = 1-1/F^2-(Rb/F/Ym)^2;
        B = 2*Rm*Rb^2/(F^2*Ym^2);
        B = B*ones(1,M);
        A = A*ones(1,M);
        con_beta = acos(sqrt(-( (Rb*Rm/F/Ym)^2.*ones(1,M) + B.^2/4./A )./R^2));
    end
else
    if F<1
        con_beta = pi/2;
    else
        A = 1-1/F^2+(Rb/F/Ym)^2;
        B = - 2*Rm*Rb^2/(F^2*Ym^2);
        B = B*ones(1,M);
        A = A*ones(1,M);
        con_beta = acos(sqrt(((Rb*Rm/F/Ym)^2.*ones(1,M) - B.^2/4./A)./R^2));
    end
end
Upper = con_beta;
end