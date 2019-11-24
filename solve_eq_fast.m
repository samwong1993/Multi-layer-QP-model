%% Zeroth-Order Algorithm for beta subprobem with estimated P/D and their gradient
%created by Huang Sen
%Email: huangsen1993@gmail.com
function beta = solve_eq_fast(R,beta,XYZ,x,k,n,p_P,p_D)
ss = 0.001;
for i = 1:100
[P gradP D gradD] = build_PD(p_P,p_D,n,beta);
dif_old = norm(XYZ(k,:)-x,2)^2 - 4*R^2*(1 - cos(D(k)/R))/2;
beta(k) = beta(k)+ss; 
[P gradP D gradD] = build_PD(p_P,p_D,n,beta);
dif = norm(XYZ(k,:)-x,2)^2 - 4*R^2*(1 - cos(D(k)/R))/2;
if sign(dif)~=sign(dif_old)
    ss = - ss/10;
else if abs(dif)>abs(dif_old)
    ss = -ss;
    end
end
if abs(dif)<1e-5
    break
end

end