%% Calculate range for flying angle
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [beta0] = rangebeta(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f)
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
beta0 = Upper+0.0001;
for iter = 1:1000
betaES = beta0;
[A B C] = QP_ABC(R,RmES,RbES,YmES,FES,betaES,0);
%Joining layer
Rmj = rm2;
Ymj = ym2;
Rbj = Rmj - Ymj;
fcj = fc2;
fj = f;
Fj = fj/fcj;
[Lower] = penetrate(Rmj,Rbj,Fj,Ymj,R,1);
betaj = Lower + 0.1;
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
Lower = acos(sqrt(-((B^2-(2*A*Rbj+B)^2)/4/A+(Rbj*Rmj/Fj/Ymj)^2)/R^2));
betaj = Lower + 0.000001;
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
[betaj tol] = angle(R,RmES,RbES,FES,YmES,betaES,0,Rmj,Rbj,Fj,Ymj,betaj,1);
%F layer
RmF = rm3;
YmF = ym3;
RbF = RmF - YmF;
fcF = fc3;
fF = f;
FF = fF/fcF;
[~,~,Lower] = beta_bound(1,FF,R,RbF,RmF,YmF);
betaF = Lower;
[betaF tol]= angle(R,RbF,Rbj,Fj,Ymj,betaj,1,RmF,RbF,FF,YmF,betaF,0);
if abs(betaF-Lower)>1e-5&tol<1e-5
    beta0 = 1.01*beta0;
end
if abs(betaF-Lower)<1e-5&tol>1e-5
    beta0 = 0.999*beta0;
end
if abs(betaF-Lower)<1e-5&tol<1e-5
    break;
end
end   
end


