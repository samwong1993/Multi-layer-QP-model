function [beta2 tol]= angle(R,Rm1,Rb1,F1,Ym1,beta_1,joining1,Rm2,Rb2,F2,Ym2,beta_2,joining2)
M = length(beta_1);
if joining1 == 0
    beta4 = pi/2*ones(1,M);
    beta1 = beta_2;
else
    beta4 = beta_2;
    beta1 = zeros(1,M);
end
acc = 0.4;
r = 0.618;
[A B C] = QP_ABC(R,Rm1,Rb1,Ym1,F1,beta_1,joining1);
delta1 = ionosphere(R,Rb1,Rm1,Ym1,A,B,C,beta_1,joining1) - ionosphere(R,Rb1,Rm1-acc,Ym1,A,B,C,beta_1,joining1);
for i = 1:5000
beta2 = r*beta1 + (1- r)*beta4;
beta3 = r*beta4 + (1- r)*beta1;
[A B C] = QP_ABC(R,Rm2,Rb2,Ym2,F2,beta1,joining2);
delta2 = ionosphere(R,Rb2,Rm1+acc,Ym2,A,B,C,beta1,joining2) - ionosphere(R,Rb2,Rm1,Ym2,A,B,C,beta1,joining2);
D1 = abs(delta2 - delta1);
[A B C] = QP_ABC(R,Rm2,Rb2,Ym2,F2,beta2,joining2);
delta2 = ionosphere(R,Rb2,Rm1+acc,Ym2,A,B,C,beta2,joining2) - ionosphere(R,Rb2,Rm1,Ym2,A,B,C,beta2,joining2);
D2 = abs(delta2 - delta1);
[A B C] = QP_ABC(R,Rm2,Rb2,Ym2,F2,beta3,joining2);
delta2 = ionosphere(R,Rb2,Rm1+acc,Ym2,A,B,C,beta3,joining2) - ionosphere(R,Rb2,Rm1,Ym2,A,B,C,beta3,joining2);
D3 = abs(delta2 - delta1);
[A B C] = QP_ABC(R,Rm2,Rb2,Ym2,F2,beta4,joining2);
delta2 = ionosphere(R,Rb2,Rm1+acc,Ym2,A,B,C,beta4,joining2) - ionosphere(R,Rb2,Rm1,Ym2,A,B,C,beta4,joining2);
D4 = abs(delta2 - delta1);
index = D2<D3;
beta4(index) = beta3(index);
dif_index = ~index;
beta1(dif_index) = beta2(dif_index);
tol = norm(delta2 - delta1);
if norm(D1-D2)<1e-5&&tol<1e-5
    break
end
end
beta = (beta1+beta2)/2;
end