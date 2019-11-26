%% Demo of Generalized Projected Gradient Descent Algorithm for Multi-layer ionospheric model
%created by Huang Sen
%Email: huangsen1993@gmail.com
clear
fid=fopen('realdata_log.txt','a+');
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
M = 5;
n = 20;
[p_P p_D] = model(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,n);
[min_dis max_dis Lower Upper] = range(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f);
%Hong Kong
[x0 y0 z0] = LGLTtoXYZ(114.16,22.28,R);
emitter = [x0 y0 z0]';
%Chong Qing
[x0 y0 z0] = LGLTtoXYZ(106.91,29.43,R);
XYZ(1,:) = [x0 y0 z0];
%Wu Han
[x0 y0 z0] = LGLTtoXYZ(114.31,30.59,R);
XYZ(2,:) = [x0 y0 z0];
%Shang Hai
[x0 y0 z0] = LGLTtoXYZ(121.47,31.23,R);
XYZ(3,:) = [x0 y0 z0];
%Xi'an
[x0 y0 z0] = LGLTtoXYZ(108.94,34.34,R);
XYZ(4,:) = [x0 y0 z0];
%Kun Ming
[x0 y0 z0] = LGLTtoXYZ(102.83,24.88,R);
XYZ(5,:) = [x0 y0 z0];
beta0 = [0.363735916327494,0.462185130143162,0.303344375146732,0.252043292085888,0.317003312557548];
% [emitter,XYZ,beta0] = generator(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,max_dis,min_dis,Lower,Upper);
beta = beta0;
x = emitter';
% [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
% for k = 1:M
%     penalty(k) = norm(XYZ(k,:)-x,2)^2 - 4*R^2*(1 - cos(D(k)/R))/2;
% end
N = M*(M-1)/2;
Omega = 0.5*(ones(N,N) + eye(N));
inv_Omega = Omega^-1;
[G] = generate_G(N,M);
sigma = [0:100:1000];
for index = 1
for noise_level = 1%1:length(sigma)
sigma_t = sigma(noise_level)* 10^-9 * 3 * 10^5 ;
noise_t0 = randn(M,1);
noise_t = (sigma_t*G*noise_t0)';
tau = generate_tau(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta0) + noise_t;
%[tau] = generate_tau(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
% [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
% obj = (G*P'-tau')'*inv_Omega*(G*P'-tau');
[x beta obj] = GPGD(M,R,G,tau,inv_Omega,Lower,Upper,max_dis,min_dis,XYZ,n,p_P,p_D);
dis = norm(x - emitter');
fprintf(fid,"%2.2f,%d,%2.2f\n",noise_level,index,dis);
end
end
pltMCCRLB('realdata_log.txt','*k-')