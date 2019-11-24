%% Function for Calculate D
%created by Huang Sen
%Email: huangsen1993@gmail.com
function output = lnD(A,B,C,r)
    output = 2*sqrt(C).*(sqrt(real(A.*r.^2 + B.*r + C))./r + B/2./sqrt(C) + sqrt(C)./r);
end