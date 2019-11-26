%% Generalized Projected Gradient Descent Algorithm for Multi-layer ionospheric model
%created by Huang Sen
%Email: huangsen1993@gmail.com
%mo is the momentum parameter and distance
function [x beta obj] = GPGD(M,R,G,tau,inv_Omega,Lower,Upper,max_dis,min_dis,XYZ,n,p_P,p_D)
        thres = 20;
        mo = 1.1;
        ss = 1e-2;
        x = [];
        x = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,R,Upper,Lower,x,n,p_P,p_D);
        beta = 0.5*(Lower + Upper)*ones(1,M);
        for i =1:5
            for k = 1:M
                 beta = solve_eq_fast(R,beta,XYZ,x,k,n,p_P,p_D);
            end
        end
        beta(beta<Lower) = Lower;
        beta(beta>Upper) = Upper;
        iter_old = 1;
        obj_min = 9999999;
        x_min = [];
        for  iter = 1:1000000
            if mod(iter,20)==0&&obj>50000
                x = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,R,Upper,Lower,x,n,p_P,p_D);
            end
            [P gradP D gradD] = build_PD(p_P,p_D,n,beta);
            obj_old = (G*P'-tau')'*inv_Omega*(G*P'-tau');
            S = XYZ;
            dBeta = zeros(M, 2);
            for i = 1:M
                dBeta(i,1) = -(S(i,1) - S(i,3)*x(1)/x(3))/R/sin(D(i)/R)/gradD(i);
                dBeta(i,2) = -(S(i,2) - S(i,3)*x(2)/x(3))/R/sin(D(i)/R)/gradD(i);
            end
            %Newton method
%             [hessP, hessD] = hessPD(A,B,C,beta,R,Rb);
%             invH = diag(1./diag(hessP));
%             dP_x = (2*G'*inv_Omega*(G*P'-tau'))'.*(invH*graP')'*dBeta;
            %Gradient descent
            dP_x = (2*G'*inv_Omega*(G*P'-tau'))'.*gradP*dBeta;
            dP_x(3) = 0;
            while(1)
                if ss*max(abs(dP_x))>thres
                    ss = 0.5*ss;
                else
                    break
                end
            end
            x = x - ss*dP_x;
            if isnan(x(1))
                break
            end
            ss = ss*mo;
            if norm(x(1:2))>=R
                x(x>R) = R;
                x(x<-R) = -R;
                x(1:2) = x(1:2)/norm(x(1:2))*R*0.99;
            end
            x(3) = sqrt(R^2 - (x(1)^2+x(2)^2));
            for i =1:5
                for k = 1:M               
                    beta = solve_eq_fast(R,beta,XYZ,x,k,n,p_P,p_D);
                end
            end
            if ~all(beta>Lower)
                if norm(x(1:2))>=R
                    x(1:2) = x(1:2)/norm(x(1:2))*R;
                end
                x(3) = - sqrt(R^2 - (x(1)^2+x(2)^2));
                for i =1:5
                    for k = 1:M
                        beta = solve_eq_fast(R,beta,XYZ,x,k,n,p_P,p_D);
                    end
                end
            end
            beta(beta<Lower) = Lower;
            beta(beta>Upper) = Upper;
            [P gradP D gradD] = build_PD(p_P,p_D,n,beta);
            obj = (G*P'-tau')'*inv_Omega*(G*P'-tau');
            if obj<1000
                thres = 5;
                mo = 1.05;
            end
            fprintf("obj:%2.2f step size:%2.6f\n",obj,ss);        
            if obj_min>obj
                iter_old = iter;
                obj_min = obj;
                x_min = x;
            end
            if abs(obj)<1e-7||(abs(iter-iter_old)>20&obj<20000)
                x = x_min;
                obj = obj_min;
                break
            end
        end
        fprintf("obj:%2.2f step size:%2.6f\n",obj,ss);
end