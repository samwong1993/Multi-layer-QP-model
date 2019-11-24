%% Calculate range for flying angle and effective detection region
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [min_dis max_dis Lower Upper] = range(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f)
RmES = rm1;
YmES = ym1;
RbES = RmES - YmES;
fcES = fc1;
fES = f;
FES = fES/fcES;
Upper = rangebeta(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f);
beta = Upper;
[P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
d = 2* R*sin(D/2/R);
min_dis = d;
[Lower] = penetrate(RmES,RbES,FES,YmES,R,0);
Lower = 1.1*Lower;
beta = Lower;
[P D] = SumPD(R,fc1,fc2,fc3,rm1,rm2,rm3,rb1,rb2,rb3,ym1,ym2,ym3,f,beta);
d = 2* R*sin(D/2/R);
max_dis = d;
end