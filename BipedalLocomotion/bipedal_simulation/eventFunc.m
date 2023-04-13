function [value, isterminal, direction] = eventFunc(t,x)
p2 = p2Func(x);
if (p2(1)> 0)
    value = p2(2);
    isterminal = 1;
    direction = -1;
else
    value = 1;
    isterminal = 0;
    direction = 0;
end

end

