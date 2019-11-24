%% Function for Calculate P
%created by Huang Sen
%Email: huangsen1993@gmail.com
function output = lnP(A,B,C,X,r)
    output = A.*r +sqrt(A.*X) +B/2;
end