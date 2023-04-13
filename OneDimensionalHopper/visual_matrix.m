function visual_matrix(T,height_sweep)

figure()
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
s=surf(H1,H2,T(2:end,2:end));
s.EdgeColor = 'none';
title('Exhaustive Method','Transition matrix excluding the absorbing state')
view(0,90)
axis tight
