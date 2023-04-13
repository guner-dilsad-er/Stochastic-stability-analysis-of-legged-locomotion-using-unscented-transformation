function [L2fh] = L2fhFunc(x,a,tp,tm)
hDotJacVel = hDotJacVelFunc(x,a,tp,tm);
Ds = DsFunc(x);
Cs = CsFunc(x);
Gs = GsFunc(x);
dhdq = dhdqFunc(x,a,tp,tm);
DsInv = inv(Ds);

L2fh = hDotJacVel + dhdq*(DsInv*(-Cs*x(6:10)-Gs));
end

