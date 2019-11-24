%% Generate tau
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [tau] = generate_tau(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta)
[P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
%Generate tau
k = 1;
for i = 1:M
    for j = i+1:M
        tau(k) = P(i) - P(j);
        k = k + 1;
    end
end
end