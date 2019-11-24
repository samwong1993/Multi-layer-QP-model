%% Function for Calculate D
%created by Huang Sen
%Email: huangsen1993@gmail.com
function output = arcsinD(A,B,C,r)
    output = asin((B.*r + 2*C)./r./sqrt( B.^2 - 4*A.*C));
end