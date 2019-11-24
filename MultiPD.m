%% Calculate P/D (E + F + Joining layer)
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [P D] = MultiPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta)
M = length(beta);
a1 = fc1^2;
a2 = a1;
a3 = fc3^2;
b1 = a1*(rb1/ym1)^2;
b3 = a3*(rb3/ym3)^2;
rc = rm3*b3*(rm3/rm1-1)/(a3 - a1 + b3*(rm3/rm1-1));
b2 = -rm3*b3*(1-rm3/rc)/rm1/(1-rm1/rc);
%Operating frequency

%Es layer
RmES = rm1;
YmES = ym1;
RbES = RmES - YmES;
fcES = fc1;
fES = f;
FES = fES/fcES;
[Upper] = penetrate(RmES,RbES,FES,YmES,R,0);
betaES = beta;
[A B C] = QP_ABC(R,RmES,RbES,YmES,FES,betaES,0);
[PES DES] = ionosphere(R,RbES,RmES,YmES,A,B,C,betaES,0);
%Joining layer
Rmj = rm2;
Ymj = ym2;
Rbj = Rmj - Ymj;
%Rj = Rm + 0.9*Ym;
fcj = fc2;
fj = f;
Fj = fj/fcj;
[Lower] = penetrate(Rmj,Rbj,Fj,Ymj,R,1);
betaj = 0.5*ones(1,M);
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
Lower = acos(sqrt(-((B.^2-(2*A*Rbj+B).^2)/4./A+(Rbj*Rmj/Fj/Ymj)^2)/R^2));
betaj = Lower + 0.000001;
% beta = [Upper+0.001:0.0001:pi/2];
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
betaj = angle(R,RmES,RbES,FES,YmES,betaES,0,Rmj,Rbj,Fj,Ymj,betaj,1);
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
[Pj Dj] = ionosphere(R,Rbj,Rmj,Ymj,A,B,C,betaj,1);
%F layer
RmF = rm3;
YmF = ym3;
RbF = RmF - YmF;
fcF = fc3;
fF = f;
FF = fF/fcF;
[Lower] = penetrate(RmF,RbF,FF,YmF,R,0);
betaF = Lower*ones(1,M);
betaF = angle(R,RbF,Rbj,Fj,Ymj,betaj,1,RmF,RbF,FF,YmF,betaF,0);
[A B C] = QP_ABC(R,RmF,RbF,YmF,FF,betaF,0);
[PF DF] = ionosphere(R,RbF,RmF,YmF,A,B,C,betaF,0);
P = PES + Pj + PF;
D = DES + Dj + DF;
% fprintf("betaES: %2.4f,betaj: %2.4f,betaF: %2.4f\n",betaES,betaj,betaF);
%save ionosphere