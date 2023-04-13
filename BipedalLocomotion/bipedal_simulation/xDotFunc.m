function [xDot] = xDotFunc(t,x,a,tp,tm)

Ds = DsFunc(x);
Cs = CsFunc(x);
Gs = GsFunc(x);
u = u96Func(x,a,tp,tm);
U=u;
if(abs(U(1)) > 75)
    U(1) = 75 * sign(U(1));
end
if(abs(U(2)) > 75)
    U(2) = 75 * sign(U(2));
end
if(abs(U(3)) > 50)
    U(3) = 50 * sign(U(3));
end
if(abs(U(4)) > 50)
    U(4) = 50 * sign(U(4));
end
u=U;
Bs = BsFunc();

xDot = [x(6:10); Ds\(-Cs*x(6:10) - Gs + Bs*u)];

end