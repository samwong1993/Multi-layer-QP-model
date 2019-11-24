%% Zeroth-Order Algorithm for beta subproblem for Multi-layer ionospheric model
%created by Huang Sen
%Email: huangsen1993@gmail.com
function beta = solve_eq(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta,XYZ,x,k,Lower,Upper)
ss = 0.001;
betaL = Lower;
betaU = Upper;
dis = norm(XYZ(k,:)-x,2);
beta0 = beta;
for i = 1:100
    [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,betaL);
    dL = 2* R*sin(D/2/R);
    [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,betaU);
    dU = 2* R*sin(D/2/R);
    betaM = 0.5*(betaL + betaU);
    [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,betaM);
    dM = 2* R*sin(D/2/R);
    if dM>dis
        betaL = betaM;
    else 
        betaU = betaM;
    end
    if abs(betaU-betaL)<1e-9
        beta = beta0;
        beta(k) = 0.5*(betaL + betaU);
        break
    end
end