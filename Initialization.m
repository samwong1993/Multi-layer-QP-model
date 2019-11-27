function x_min = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,R,Upper,Lower,x_input,n,p_P,p_D)  
while(1)
    x = randn(1,3);
    x = x/norm(x)*R;
    dis = [];
    for iter = 1:M
        dis = [dis norm(XYZ(iter,:) - x)];
    end
    if all(dis<max_dis)&all(dis>min_dis)
        break
    end
end
x0 = x;
beta = 0.5*(Lower + Upper)*ones(1,M);
for i =1:5
    for k = 1:M
        beta = solve_eq_fast(R,beta,XYZ,x,k,n,p_P,p_D,Lower,Upper);
    end
end
beta(beta<Lower) = Lower;
beta(beta>Upper) = Upper;
[P gradP D gradD] = build_PD(p_P,p_D,n,beta);
obj_x0 = (G*P'-tau')'*inv_Omega*(G*P'-tau');
center = R*sum(XYZ)/M/(norm(sum(XYZ)/M));
x_ini = 2*(center*(x'/R)*center/R - x) + x;
for i =1:5
    for k = 1:M
        beta = solve_eq_fast(R,beta,XYZ,x_ini,k,n,p_P,p_D,Lower,Upper);
    end
end
beta(beta<Lower) = Lower;
beta(beta>Upper) = Upper;
[P gradP D gradD] = build_PD(p_P,p_D,n,beta);
obj_ini = (G*P'-tau')'*inv_Omega*(G*P'-tau');
if obj_ini>obj_x0
    obj_min = obj_x0;
    x_min = x0;
else
    obj_min = obj_ini;
    x_min = x_ini;
end

if ~isempty(x_input)
    beta = 0.5*(Lower + Upper)*ones(1,M);
    for i =1:5
        for k = 1:M
            beta = solve_eq_fast(R,beta,XYZ,x_input,k,n,p_P,p_D,Lower,Upper);
        end
    end
    beta(beta<Lower) = Lower;
    beta(beta>Upper) = Upper;
    [P gradP D gradD] = build_PD(p_P,p_D,n,beta);
    obj_input = (G*P'-tau')'*inv_Omega*(G*P'-tau');
    if obj_min>obj_input
        obj_min = obj_input;
        x_min = x_input;
    end
end

end