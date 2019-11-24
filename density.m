clear
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
r = 6465:1:6650;
for i = 1 :length(r)
    if r(i) <= 6485
        obj(i) = a1 - b1*(1 - rm1/r(i))^2;
    else if r(i) <= rc
            obj(i) = a2 + b2*(1 - rm1/r(i))^2;
        else if r(i) <= 6650
                obj(i) = a3 - b3*(1 - rm3/r(i))^2;
            end
        end
    end
end
plot(obj,r)