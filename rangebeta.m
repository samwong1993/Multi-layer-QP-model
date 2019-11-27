%% Calculate range for flying angle
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [betaL] = rangebeta(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f)
a1 = fc1^2;
a2 = a1;
a3 = fc3^2;
b1 = a1*(rb1/ym1)^2;
b3 = a3*(rb3/ym3)^2;
rc = rm3*b3*(rm3/rm1-1)/(a3 - a1 + b3*(rm3/rm1-1));
b2 = -rm3*b3*(1-rm3/rc)/rm1/(1-rm1/rc);
%Es layer
RmES = rm1;
YmES = ym1;
RbES = RmES - YmES;
fcES = fc1;
fES = f;
FES = fES/fcES;
[Upper] = penetrate(RmES,RbES,FES,YmES,R,0);
betaL = Upper+0.0001;
betaU = pi/2;
betaM = 0.5*(betaL + betaU);
for iter = 1:1000
betaM = 0.5*(betaL + betaU);
[betaFM tolM LowerM] = bisection(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,betaM);
if abs(betaFM-LowerM)>1e-5&tolM<1e-5
    betaL = betaM;
end
if abs(betaFM-LowerM)<1e-5&tolM>1e-5
    betaU = betaM;
end
if abs(betaFM-LowerM)<1e-5&tolM<1e-5
    break;
end
end   
end


