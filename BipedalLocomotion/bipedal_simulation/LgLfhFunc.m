function [LgLfh] = LgLfhFunc(x,a,tp,tm)
dhdq = dhdqFunc(x,a,tp,tm);
Ds = DsFunc(x);
DsInv = inv(Ds);
Bs = BsFunc();

LgLfh = dhdq*DsInv*Bs;
end

