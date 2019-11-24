%% Calculate estimated P/D and their gradient (E + F + Joining + Free space)
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [P gradP D gradD] = build_PD(p_P,p_D,n,beta)
    M = length(beta);
    P = polyval(p_P,beta);
    pp = polyder(p_P);
    gradP = polyval(pp,beta);
% 	P = zeros(1,M);
%     for j = 1:n+1
%         P = P + p_P(j)*beta.^(n+1-j);
%     end
%     gradP = zeros(1,M);
% 	for j = 1:n
%         gradP = gradP + p_P(j)*(n+1-j)*beta.^(n-j);
%     end
    D = polyval(p_D,beta);
	dd = polyder(p_D);
    gradD = polyval(dd,beta);
% 	D = zeros(1,M);
%     for j = 1:n+1
%         D = D + p_D(j)*beta.^(n+1-j);
%     end
%     gradD = zeros(1,M);
% 	for j = 1:n
%         gradD = gradD + p_D(j)*(n+1-j)*beta.^(n-j);
%     end
end