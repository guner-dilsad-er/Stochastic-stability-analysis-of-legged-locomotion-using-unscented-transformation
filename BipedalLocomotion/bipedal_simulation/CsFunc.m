function Cs = CsFunc(in1)
%CSFUNC
%    CS = CSFUNC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    18-Dec-2021 16:49:21

dq1 = in1(6,:);
dq2 = in1(7,:);
dq3 = in1(8,:);
dq4 = in1(9,:);
dq5 = in1(10,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
Cs = reshape([dq3.*sin(q3).*(-4.3088)+5.0,dq1.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))+dq5.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))+dq3.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)),sin(q3).*(dq1+dq5).*4.3088,dq3.*sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+dq1.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))+dq5.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)),dq3.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)+sin(q1+q3).*(1.44e+2./1.25e+2)-sin(q3).*4.3088)+dq1.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+sin(q1+q3).*(1.44e+2./1.25e+2)+sin(q1).*(1.44e+2./1.25e+2))+dq5.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+sin(q1+q3).*(1.44e+2./1.25e+2)+sin(q1).*(1.44e+2./1.25e+2)),-dq2.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))-dq5.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))-dq4.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)),dq4.*sin(q4).*(-1.92e+2./6.25e+2)+5.0,dq4.*sin(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-dq2.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2))-dq5.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)),sin(q4).*(dq2+dq5).*(1.92e+2./6.25e+2),-dq4.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+sin(q4).*(1.92e+2./6.25e+2))-dq2.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))-dq5.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)),sin(q3).*(dq1+dq3+dq5).*(-4.3088),(sin(q1-q2+q3-q4).*6.4e+1+sin(q1-q2+q3).*1.69e+2).*(dq1+dq3+dq5).*(3.0./6.25e+2),1.0e+1,sin(q1-q2+q3-q4).*(dq1+dq3+dq5).*(1.92e+2./6.25e+2),((dq1+dq3+dq5).*(sin(q1-q2+q3-q4).*1.92e+2+sin(q1-q2+q3).*5.07e+2+sin(q1+q3).*7.2e+2-sin(q3).*2.693e+3))./6.25e+2,(sin(q1-q2+q3-q4)-sin(-q1+q2+q4)).*(dq2+dq4+dq5).*(-1.92e+2./6.25e+2),sin(q4).*(dq2+dq4+dq5).*(-1.92e+2./6.25e+2),sin(q1-q2+q3-q4).*(dq2+dq4+dq5).*(-1.92e+2./6.25e+2),1.0e+1,(sin(q1-q2+q3-q4)-sin(-q1+q2+q4)+sin(q4)).*(dq2+dq4+dq5).*(-1.92e+2./6.25e+2),dq2.*sin(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-dq4.*sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-dq5.*sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-dq2.*sin(q1-q2+q3).*(5.07e+2./6.25e+2)+dq2.*sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+dq4.*sin(-q1+q2+q4).*(1.92e+2./6.25e+2)-dq5.*sin(q1-q2+q3).*(5.07e+2./6.25e+2)+dq5.*sin(-q1+q2+q4).*(1.92e+2./6.25e+2)-dq5.*sin(q1+q3).*(1.44e+2./1.25e+2)-dq3.*sin(q3).*4.3088-dq5.*sin(q1).*(1.44e+2./1.25e+2)-dq2.*sin(q1-q2).*(5.07e+2./6.25e+2)-dq5.*sin(q1-q2).*(5.07e+2./6.25e+2),dq1.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))+dq5.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))-dq4.*sin(q4).*(1.92e+2./6.25e+2)+dq3.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)),-dq5.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)+sin(q1+q3).*(1.44e+2./1.25e+2)-sin(q3).*4.3088)-dq4.*sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+dq1.*sin(q3).*4.3088-dq2.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)),dq5.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+sin(q4).*(1.92e+2./6.25e+2))+dq3.*sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+dq2.*sin(q4).*(1.92e+2./6.25e+2)+dq1.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)),-dq4.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+sin(q4).*(1.92e+2./6.25e+2))+dq3.*(sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)+sin(q1+q3).*(1.44e+2./1.25e+2)-sin(q3).*4.3088)-dq2.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2))+dq1.*(sin(q1-q2).*(5.07e+2./6.25e+2)+sin(q1-q2+q3-q4).*(1.92e+2./6.25e+2)+sin(q1-q2+q3).*(5.07e+2./6.25e+2)-sin(-q1+q2+q4).*(1.92e+2./6.25e+2)+sin(q1+q3).*(1.44e+2./1.25e+2)+sin(q1).*(1.44e+2./1.25e+2))],[5,5]);
