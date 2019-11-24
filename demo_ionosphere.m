%% Demo for inonsphere model
%created by Huang Sen
%Email: huangsen1993@gmail.com
clear
R = 6371.2;
%Three layers
fc1 = 4;
fc2 = 4;
fc3 = 10;
rm1 = 6485;
rm2 = 6550;
rm3 = 6650;
rb1 = 6465;
rb2 = rm1;
rb3 = 6550;
ym1 = 20;
ym2 = 65;
ym3 = 100;
a1 = fc1^2;
a2 = a1;
a3 = fc3^2;
b1 = a1*(rb1/ym1)^2;
b3 = a3*(rb3/ym3)^2;
rc = rm3*b3*(rm3/rm1-1)/(a3 - a1 + b3*(rm3/rm1-1));
b2 = -rm3*b3*(1-rm3/rc)/rm1/(1-rm1/rc);
%Operating frequency
f = 15;
%Es layer
RmES = 6485;
YmES = 20;
RbES = RmES - YmES;
%Rj = Rm + 0.9*Ym;
fcES = fc1;
fES = f;
FES = fES/fcES;
[Upper] = penetrate(RmES,RbES,FES,YmES,R,0);
betaES = 0.4;
[A B C] = QP_ABC(R,RmES,RbES,YmES,FES,betaES,0);
for i = 6465:6485
    % for i = 6465:1:6485
    %     [P(i) D(i)] = ionosphere(R,RbES,i,YmES,A,B,C,beta,joining);
    % end
    % % Rt= -(B+sqrt(B.^2 - 4*A.*C))/2./A;
    % plot([6465:1:6485],P(6465:1:6485))
    [PES(i-6465+1) DES(i-6465+1)] = ionosphere(R,RbES,i,YmES,A,B,C,betaES,0);
end

%Joining layer
Rmj = rm2;
Ymj = ym2;
Rbj = Rmj - Ymj;
%Rj = Rm + 0.9*Ym;
fcj = fc2;
fj = f;
Fj = fj/fcj;
[Lower] = penetrate(Rmj,Rbj,Fj,Ymj,R,1);
betaj = 0.5;
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
Lower = acos(sqrt(-((B^2-(2*A*Rbj+B)^2)/4/A+(Rbj*Rmj/Fj/Ymj)^2)/R^2));
betaj = Lower + 0.000001;
% beta = [Upper+0.001:0.0001:pi/2];
[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
% gradP(A,B,C,R,betaj,6485,Rbj,joining)
% for i = 6485:1:6550
% [P(i) D(i)] = ionosphere(R,Rbj,i,Ymj,A,B,C,betaj,joining);
% end
% % Rt= -(B+sqrt(B.^2 - 4*A.*C))/2./A;
% plot([6485:1:6550],P(6485:1:6550))
betaj = angle(R,RmES,RbES,FES,YmES,betaES,0,Rmj,Rbj,Fj,Ymj,betaj,1);

[A B C] = QP_ABC(R,Rmj,Rbj,Ymj,Fj,betaj,1);
for i = 6486:6550
    [Pj(i-6486+1) Dj(i-6486+1)] = ionosphere(R,Rbj,i,Ymj,A,B,C,betaj,1);
end

%F layer
RmF = rm3;
YmF = ym3;
RbF = RmF - YmF;
fcF = fc3;
fF = f;
FF = fF/fcF;
[Lower] = penetrate(RmF,RbF,FF,YmF,R,0);
betaF = Lower;
betaF = angle(R,RbF,Rbj,Fj,Ymj,betaj,1,RmF,RbF,FF,YmF,betaF,0);
[A B C] = QP_ABC(R,RmF,RbF,YmF,FF,betaF,0);
Rt = -(B + sqrt(B^2-4*A*C))/2/A;
for i = 6551:floor(Rt)
    Rt = i;
    DF(i-6551+1) = - R^2*cos(betaF)./sqrt(C).*log(abs(lnD(A,B,C,Rt))./lnD(A,B,C,RbF));
    Xb = abs(A.*RbF.^2 + B.*RbF + C);
    Xm = abs(A.*Rt.^2 + B.*Rt + C);
    PF(i-6551+1) = ( sqrt(Xm) - sqrt(Xb) )/A + B/2/sqrt(A^3)*log(lnP(A,B,C,Xb,RbF)/lnP(A,B,C,Xm,Rt));
end
P = [PES max(PES)+Pj max(PES)+max(Pj) + PF];
D = [DES max(DES)+Dj max(DES)+max(Dj) + DF];
plot(P,[6465:floor(Rt)], 'g-', 'linewidth', 1.5)
title('Demo of Three Layers Ionosphere Model')
xlabel('Distance')
ylabel('Height')
str = ['\beta = ' num2str(betaES)];
legend(str)
hold on
fprintf("%2.4f,%2.4f,%2.4f\n",betaES,betaj,betaF)