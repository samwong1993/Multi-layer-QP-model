function grad = gradP(A,B,C,R,beta,r,Rb,joining)
if joining == true
    dC_dB = - R^2*2*cos(beta).*sin(beta);
    Xb = A*Rb^2 + B*Rb + C;
    dsqrt_Xb = 0.5*Xb.^-0.5.*dC_dB;
    X = A*r^2 + B*r + C;
    dsqrt_X = 0.5*X.^-0.5.*dC_dB;
    delta = B.^2 - 4*A.*C;
    tau = (2*A*r + B)./-delta;
    deltab = B.^2 - 4*A.*C;
    taub = (2*A*Rb + B)./-deltab;
    grad = 1./A.*(dsqrt_X - dsqrt_Xb) + B/2./(-A).^(1.5).*(1./sqrt(1 - tau.^2)*-2.*(A*r+B)*-0.5.*(delta).^(-1.5)*-4.*A.*dC_dB  - 1./sqrt(1 - taub.^2)*-2.*(A*Rb+B)*-0.5.*(deltab).^(-1.5)*-4.*A.*dC_dB   );
else
%     M = length(beta);
%     Rt = ones(1,M)*r;
%     index = (B.^2 - 4*A.*C>=0);
%     Rt(index) = -(B(index)+sqrt(B(index).^2 - 4*A(index).*C(index)))/2./A(index);
%     if r<Rt
%         Rt = r;
%     end
    dC_dB = R^2*2*cos(beta).*sin(beta);
    Xb = A*Rb^2 + B*Rb + C;
    dsqrt_Xb = 0.5*Xb.^-0.5.*dC_dB;
    X = abs(A.*r.^2 + B.*r + C);
    dsqrt_X = 0.5*X.^-0.5.*dC_dB;
    grad = 1./A.*(dsqrt_X - dsqrt_Xb) + B/2./A.^(1.5).*( 1./(sqrt(A.*Xb)+A*Rb+B/2).*sqrt(A).*dsqrt_Xb - 1./(sqrt(A.*X)+A.*r+B/2).*sqrt(A).*dsqrt_X);
end
end