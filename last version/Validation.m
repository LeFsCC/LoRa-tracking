syms x;
a1 = [1.1:0.1:20.1];
a2 = [];

for i = 1:length(a1)
   eq = x^2-x*(a1(i)+1)-1;
   a2 = [a2 double(solve(eq, x > 0))]; 
end

c = 1;
s = [];

for i = 1:length(a1)
    syms u v;
    b1 = sqrt(a1(i)^2-c^2);
    b2 = sqrt(a2(i)^2-c^2);
    eq1 = (u-c)^2/a2(i)^2+v^2/b2^2-1;
    eq2 = u^2/b1^2+(v+c)^2/a1(i)^2-1;
    s = [s solve(eq1,eq2,u,v)];
end

cach_u = [];
cach_v = [];
sel = [];

for i = 1:4
    if isreal(double(s(1).u(i)))
        sel = [sel i];
    end
end

for i = 1:length(s)
    for j = 1:4
        if isreal(double(s(i).u(j)))
            cach_u = [cach_u double(s(i).u(j))];
            cach_v = [cach_v double(s(i).v(j))];
        end
    end
end

for i = 1:length(cach_u)
   scatter(cach_u(1:i), cach_v(1:i));
   set(0,'defaultfigurecolor','w');
   pause(0.04);
end
