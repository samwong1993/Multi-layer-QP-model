%% Function for Calculate D
%created by Huang Sen
%Email: huangsen1993@gmail.com
function output = arcsinP(A,B,C,r)
    output = asin((2*A.*r+B)./-sqrt(B.^2-4*A.*C));
end
