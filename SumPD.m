%% Calculate P/D (E + F + Joining + Free space)
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta)
RmES = rm1;
YmES = ym1;
RbES = RmES - YmES;
fcES = fc1;
fES = f;
FES = fES/fcES;
[P D] = MultiPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
Gamma = acos(R/RbES.*cos(beta));
D0 = R*(Gamma - beta);
P0 = RbES*sin(Gamma) - R*sin(beta);
P = 2*(P + P0);
D = 2*(D + D0);
end
