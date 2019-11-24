%% Generate the ionosphere model and calculate the range for flying angle
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [p_P p_D] = model(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,n)
%Es layer
RmES = rm1;
YmES = ym1;
RbES = RmES - YmES;
fcES = fc1;
fES = f;
FES = fES/fcES;
[Upper] = penetrate(RmES,RbES,FES,YmES,R,0);
beta = [0:0.001:Upper-0.001];
for iter = 1:length(beta)
    [A B C] = QP_ABC(R,RmES,RbES,YmES,FES,beta(iter),0);
    [P(iter) D(iter)] = PD(A,B,C,beta(iter),R,RbES);
end

beta0 = [Upper+0.001:0.001:pi/2];
for iter = 1:length(beta0)
betaES = beta0(iter);
[A B C] = QP_ABC(R,RmES,RbES,YmES,FES,betaES,0);
[PES DES] = ionosphere(R,RbES,RmES,YmES,A,B,C,betaES,0);
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
[betaj, tol] = angle(R,RmES,RbES,FES,YmES,betaES,0,Rmj,Rbj,Fj,Ymj,betaj,1);
if tol>1e-5
    break
end
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
betaF = Lower;
[betaF tol]= angle(R,RbF,Rbj,Fj,Ymj,betaj,1,RmF,RbF,FF,YmF,betaF,0);
if abs(betaF - Lower)<1e-6
    break
end
[A B C] = QP_ABC(R,RmF,RbF,YmF,FF,betaF,0);
[PF DF] = ionosphere(R,RbF,RmF,YmF,A,B,C,betaF,0);
Gamma = acos(R/RbES*cos(betaES));
D0 = R*(Gamma - betaES);
P0 = RbES*sin(Gamma) - R*sin(betaES);
MultiP(iter) = 2*(P0 + PES + Pj + PF);
MultiD(iter) = 2*(D0 + DES + Dj + DF);
end
P = [P MultiP];
D = [D MultiD];
beta = [beta beta0];
figure(1)
plot(beta(1:length(P)),P, 'r--', 'linewidth', 1.5)
hold on
plot(beta(1:length(D)),D, 'b:', 'linewidth', 1.5)
xlabel('\beta')
ylabel('Function Value')
legend('P','D')
xlim([min(beta(1:length(P))) max(beta(1:length(P)))])
%fit the function P
[min_dis max_dis Lower Upper] = range(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f);
betatmp = beta(1:length(D));
index = find(betatmp<Upper&betatmp>Lower);
betatmp = betatmp(index);
Ptmp = P(index);
Dtmp = D(index);
p_P = polyfit(betatmp,Ptmp,n);
p_D = polyfit(betatmp,Dtmp,n);
for i =1:length(Ptmp)
    Pval(i) = 0;
    for j = 1:n+1
        Pval(i) = Pval(i) + p_P(j)*betatmp(i)^(n+1-j);
    end
end
for i =1:length(Dtmp)
    Dval(i) = 0;
    for j = 1:n+1
        Dval(i) = Dval(i) + p_D(j)*betatmp(i)^(n+1-j);
    end
end
figure(2)
plot(betatmp,Dtmp,'k--', 'linewidth',2)
hold on
plot(betatmp,Dval,'r-','linewidth', 0.5)
xlabel('\beta')
ylabel('Function Value')
legend('True D','Estmated D')
%save('fit.mat','betatmp','Ptmp','Dtmp')
norm(Ptmp - Pval)
norm(Dtmp - Dval)
end