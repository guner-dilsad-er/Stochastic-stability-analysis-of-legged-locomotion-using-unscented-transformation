function E2 = E2Func(in1)
%E2FUNC
%    E2 = E2FUNC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    18-Dec-2021 16:49:22

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
E2 = reshape([cos(q1+q3+q5).*(2.0./5.0)+cos(q1+q5).*(2.0./5.0),sin(q1+q3+q5).*(2.0./5.0)+sin(q1+q5).*(2.0./5.0),cos(q2+q4+q5).*(-2.0./5.0)-cos(q2+q5).*(2.0./5.0),sin(q2+q4+q5).*(-2.0./5.0)-sin(q2+q5).*(2.0./5.0),cos(q1+q3+q5).*(2.0./5.0),sin(q1+q3+q5).*(2.0./5.0),cos(q2+q4+q5).*(-2.0./5.0),sin(q2+q4+q5).*(-2.0./5.0),cos(q1+q3+q5).*(2.0./5.0)-cos(q2+q4+q5).*(2.0./5.0)+cos(q1+q5).*(2.0./5.0)-cos(q2+q5).*(2.0./5.0),sin(q1+q3+q5).*(2.0./5.0)-sin(q2+q4+q5).*(2.0./5.0)+sin(q1+q5).*(2.0./5.0)-sin(q2+q5).*(2.0./5.0),1.0,0.0,0.0,1.0],[2,7]);