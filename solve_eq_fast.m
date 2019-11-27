%% Zeroth-Order Algorithm for beta subprobem with estimated P/D and their gradient
%created by Huang Sen
%Email: huangsen1993@gmail.com
function beta = solve_eq_fast(R,beta,XYZ,x,k,n,p_P,p_D,Lower,Upper)
betaL = Lower;
betaU = Upper;
dis = norm(XYZ(k,:)-x,2);
beta0 = beta;
for i = 1:1000
    [P gradP D gradD] = build_PD(p_P,p_D,n,betaL);
    dL = 2* R*sin(D/2/R);
    [P gradP D gradD] = build_PD(p_P,p_D,n,betaU);
    dU = 2* R*sin(D/2/R);
    betaM = 0.5*(betaL + betaU);
    [P gradP D gradD] = build_PD(p_P,p_D,n,betaM);
    dM = 2* R*sin(D/2/R);
    if dM>dis
        betaL = betaM;
    else 
        betaU = betaM;
    end
    if abs(betaU-betaL)<1e-12
        beta = beta0;
        beta(k) = 0.5*(betaL + betaU);
        break
    end
end

end