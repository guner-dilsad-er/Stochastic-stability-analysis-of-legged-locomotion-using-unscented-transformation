function J=jcbn_finite_diff(f,x0)
% Adapted from Mathews, J. H., & Fink, K. D. (2004). Numerical methods using MATLAB
% f: function
% x0: the point where jacobian is evaluated

x0 = x0(:);

%h=10.^(-(6:0.5:12));

h=10.^(-(1:0.5:12));
ni = length(x0);
no = length(f(x0'));

nd = length(h);

J = zeros(no,ni);

for i = 1:ni
    index = [];
    D = zeros(no,nd);
    D_abs = zeros(1,nd);
    dummy = zeros(ni,1);
    dummy(i) = 1; % to select desired direction
    for j = 1:nd
        % D(:,j) = (-f((x0+2*h(j)*dummy)')+8*f((x0+h(j)*dummy)')-8*f((x0-h(j)*dummy)')+f((x0-2*h(j)*dummy)'))/12/h(j);
        %  D(:,j) = (-f((x0+x0*2*h(j)*dummy)')+8*f((x0+x0*h(j)*dummy)')-8*f((x0-x0*h(j)*dummy)')+f((x0-x0*2*h(j)*dummy)'))/12/h(j);
        D(:,j) = (f((x0+h(j)*dummy)')-f((x0-h(j)*dummy)'))/2/h(j);
        D_abs(j) = norm(D(:,j));
        if j > 2
            dummy1 = diff(abs(diff(D_abs(1:j))))>0;
            index = find(dummy1); % DEBUG here and see h(j)*dummy is close to edge.
            if ~isempty(index)
                break;
            end
        end
    end
    %     dummy1 = diff(abs(diff(D_abs)))>0; %% h(1) and h(end) may be the answer !!!!
    if isempty(index)
        j = 2;
    end
    J(:,i) = D(:,j-1);
end