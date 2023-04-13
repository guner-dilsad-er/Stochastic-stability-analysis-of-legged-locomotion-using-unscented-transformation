function De = DeFunc(in1)
%DEFUNC
%    DE = DEFUNC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    18-Dec-2021 16:49:21

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
De = reshape([cos(q3).*8.6176+9.4518,cos(q1-q2).*(-5.07e+2./6.25e+2)-cos(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2),cos(q3).*4.3088+4.88992,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2),cos(q1-q2).*(-5.07e+2./6.25e+2)-cos(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2)-cos(q1+q3).*(1.44e+2./1.25e+2)-cos(q1).*(1.44e+2./1.25e+2)+cos(q3).*8.6176+9.4518,cos(q1+q3+q5).*1.2032e+1+cos(q1+q5).*1.0772e+1,sin(q1+q3+q5).*1.2032e+1+sin(q1+q5).*1.0772e+1,cos(q1-q2).*(-5.07e+2./6.25e+2)-cos(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2),cos(q4).*(3.84e+2./6.25e+2)+1.4486,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2),cos(q4).*(1.92e+2./6.25e+2)+3.8432e-1,cos(q1-q2).*(-5.07e+2./6.25e+2)-cos(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2)+cos(q4).*(3.84e+2./6.25e+2)+1.4486,cos(q2+q4+q5).*(-9.6e+1./1.25e+2)-cos(q2+q5).*(5.07e+2./2.5e+2),sin(q2+q4+q5).*(-9.6e+1./1.25e+2)-sin(q2+q5).*(5.07e+2./2.5e+2),cos(q3).*4.3088+4.88992,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2),4.88992,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2),cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(q1+q3).*(1.44e+2./1.25e+2)+cos(q3).*4.3088+4.88992,cos(q1+q3+q5).*1.2032e+1,sin(q1+q3+q5).*1.2032e+1,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2),cos(q4).*(1.92e+2./6.25e+2)+3.8432e-1,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2),3.8432e-1,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2)+cos(q4).*(1.92e+2./6.25e+2)+3.8432e-1,cos(q2+q4+q5).*(-9.6e+1./1.25e+2),sin(q2+q4+q5).*(-9.6e+1./1.25e+2),cos(q1-q2).*(-5.07e+2./6.25e+2)-cos(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2)-cos(q1+q3).*(1.44e+2./1.25e+2)-cos(q1).*(1.44e+2./1.25e+2)+cos(q3).*8.6176+9.4518,cos(q1-q2).*(-5.07e+2./6.25e+2)-cos(q1-q2+q3-q4).*(1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2)+cos(q4).*(3.84e+2./6.25e+2)+1.4486,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(q1-q2+q3).*(5.07e+2./6.25e+2)-cos(q1+q3).*(1.44e+2./1.25e+2)+cos(q3).*4.3088+4.88992,cos(q1-q2+q3-q4).*(-1.92e+2./6.25e+2)-cos(-q1+q2+q4).*(1.92e+2./6.25e+2)+cos(q4).*(1.92e+2./6.25e+2)+3.8432e-1,cos(q1-q2).*(-1.6224)-cos(q1-q2+q3-q4).*(3.84e+2./6.25e+2)-cos(q1-q2+q3).*1.6224-cos(-q1+q2+q4).*(3.84e+2./6.25e+2)-cos(q1+q3).*(2.88e+2./1.25e+2)-cos(q1).*(2.88e+2./1.25e+2)+cos(q3).*8.6176+cos(q4).*(3.84e+2./6.25e+2)+1.29216e+1,cos(q1+q3+q5).*1.2032e+1-cos(q2+q4+q5).*(9.6e+1./1.25e+2)+cos(q1+q5).*1.0772e+1-cos(q2+q5).*(5.07e+2./2.5e+2)-cos(q5).*(7.2e+1./2.5e+1),sin(q1+q3+q5).*1.2032e+1-sin(q2+q4+q5).*(9.6e+1./1.25e+2)+sin(q1+q5).*1.0772e+1-sin(q2+q5).*(5.07e+2./2.5e+2)-sin(q5).*(7.2e+1./2.5e+1),cos(q1+q3+q5).*1.2032e+1+cos(q1+q5).*1.0772e+1,cos(q2+q4+q5).*(-9.6e+1./1.25e+2)-cos(q2+q5).*(5.07e+2./2.5e+2),cos(q1+q3+q5).*1.2032e+1,cos(q2+q4+q5).*(-9.6e+1./1.25e+2),cos(q1+q3+q5).*1.2032e+1-cos(q2+q4+q5).*(9.6e+1./1.25e+2)+cos(q1+q5).*1.0772e+1-cos(q2+q5).*(5.07e+2./2.5e+2)-cos(q5).*(7.2e+1./2.5e+1),3.2e+1,0.0,sin(q1+q3+q5).*1.2032e+1+sin(q1+q5).*1.0772e+1,sin(q2+q4+q5).*(-9.6e+1./1.25e+2)-sin(q2+q5).*(5.07e+2./2.5e+2),sin(q1+q3+q5).*1.2032e+1,sin(q2+q4+q5).*(-9.6e+1./1.25e+2),sin(q1+q3+q5).*1.2032e+1-sin(q2+q4+q5).*(9.6e+1./1.25e+2)+sin(q1+q5).*1.0772e+1-sin(q2+q5).*(5.07e+2./2.5e+2)-sin(q5).*(7.2e+1./2.5e+1),0.0,3.2e+1],[7,7]);
