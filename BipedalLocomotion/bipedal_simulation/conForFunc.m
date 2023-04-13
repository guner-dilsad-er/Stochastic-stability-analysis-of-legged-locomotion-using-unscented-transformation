function [F] = conForFunc(x,a,tp,tm)

De = DeFunc(x);
Ce = CeFunc(x);
Ge = GeFunc(x);

xDot = xDotFunc(0,x,a,tp,tm);
qDot = [xDot(1:5);zeros(2,1)];
qDDot = [xDot(6:10);zeros(2,1)];

F = De(end-1:end,:)*qDDot + Ce(end-1:end,:)*qDot + Ge(end-1:end);
end

