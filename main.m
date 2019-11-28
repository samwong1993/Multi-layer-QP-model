%% Demo of Generalized Projected Gradient Descent Algorithm for Multi-layer ionospheric model
%created by Huang Sen
%Email: huangsen1993@gmail.com
clear
fid=fopen('realdata_log.txt','a+');
R = 6371.2;
%plt = 1 plot the earth,emitter,sensors and the generated sequence
plt = 0;
if plt == 1
    figure('color','k')
    %Add legend
    %emitter
    point1 = scatter3(0,0,0,50,'filled','r');
    hold on
    %sensors
    point2 = scatter3(0,0,0,'filled','c');
    %initial point
    point3 = scatter3(0,0,0,50,'filled','b');
    %recover location
    point4 = scatter3(0,0,0,50,'k*');
    %generated sequence
    point5 = scatter3(0,0,0,5,'m');
    %plot earth
    npanels = 72;
    alpha   = 1;
    image_file = 'earth.jpg';
    [x0, y0, z0] = ellipsoid(0, 0, 0, R, R, R);
    globe = surf(x0, y0, -z0, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    cdata = imread(image_file);
    set(gca, 'NextPlot','add', 'Visible','off');
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    axis off;
    axis equal;
    axis auto;
end
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
beta0 = [0.363735916581471,0.462185130510566,0.303344375159694,0.252043292174790,0.317003313068478];
%[emitter,XYZ,beta0] = generator(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,max_dis,min_dis,Lower,Upper);
beta = beta0;
x = emitter';
% [P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
% for k = 1:M
%     penalty(k) = norm(XYZ(k,:)-x,2) - 2* R*sin(D(k)/2/R);
% end
% penalty
% obj = (G*P'-tau')'*inv_Omega*(G*P'-tau');
N = M*(M-1)/2;
Omega = 0.5*(ones(N,N) + eye(N));
inv_Omega = Omega^-1;
[G] = generate_G(N,M);
tau = generate_tau(M,R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta0);
[x beta obj] = GPGD(M,R,G,tau,inv_Omega,Lower,Upper,max_dis,min_dis,XYZ,n,p_P,p_D,plt);
dis = norm(x - emitter');
fprintf("The distance to the emitter is: %2.2f km\n",dis)
fprintf("%s\n",'Localization successful!')
fprintf("True Location:(%2.2f,%2.2f,%2.2f)\n",emitter(1),emitter(2),emitter(3))
fprintf("Location:(%2.2f,%2.2f,%2.2f)\n",x(1),x(2),x(3))
%Output results
if plt == 1
    scatter3(emitter(1),emitter(2),emitter(3),50,'filled','r')
    scatter3(x(1),x(2),x(3),50,'k*')
    %text(emitter(1),emitter(2),emitter(3),'e')
end
for i = 1:M
    if plt == 1
        scatter3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'filled','c')
        temp = strcat('s ',num2str(i));
        text(XYZ(i,1),XYZ(i,2),XYZ(i,3),temp);
    end
    fprintf("Sensor %d:(%2.2f,%2.2f,%2.2f)\n",i,XYZ(i,1),XYZ(i,2),XYZ(i,3))
end
fprintf("True Flying angle:(")
for i = 1:M-1
    fprintf("%2.2f,",beta0(i))
end
fprintf("%2.2f)\n",beta0(M))
fprintf("Flying angle:(")
for i = 1:M-1
    fprintf("%2.2f,",beta(i))
end
fprintf("%2.2f)\n",beta(M))
%Lengend
if plt == 1
	h = legend([point1,point2,point3,point4,point5],'Emitter', 'Sensors', 'Initial Point','Estimated Location','Generated Points','AutoUpdate','off');
    %set(h,'box','off')
end