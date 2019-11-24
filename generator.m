%% Generate emitter and sensors for Multi-layer ionospheric model
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [emitter,XYZ,beta0] = generator(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,max_dis,min_dis,Lower,Upper)
emitter = randn(3,1);
emitter = R*emitter / norm(emitter);
for i = 1:M
while(1)
    XYZ(i,:) = randn(1,3);
    XYZ(i,:) = R*XYZ(i,:) / norm(XYZ(i,:));
    if norm(emitter' - XYZ(i,:))<max_dis&norm(emitter' - XYZ(i,:))>min_dis
        break
    end
end
end
%Calculate corresponing beta of x
while(1)
beta = 0.5*(Lower+Upper)*ones(1,M);
x = emitter';
for i =1:5
    for k = 1:M
        beta = solve_eq(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta,XYZ,x,k,Lower,Upper);
    end
end
[P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
for k = 1:M
    penalty(k) = norm(XYZ(k,:)-x,2) - 2* R*sin(D(k)/2/R);
end
beta0 = beta;
if sum(penalty) < 1e-2
    break;
end
end
end