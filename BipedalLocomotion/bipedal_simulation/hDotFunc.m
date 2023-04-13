function hDot = hDotFunc(in1,in2,tp,tm)
%HDOTFUNC
%    HDOT = HDOTFUNC(IN1,IN2,TP,TM)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    18-Dec-2021 16:49:23

a1_1 = in2(1);
a1_2 = in2(5);
a1_3 = in2(9);
a1_4 = in2(13);
a1_5 = in2(17);
a1_6 = in2(21);
a1_7 = in2(25);
a2_1 = in2(2);
a2_2 = in2(6);
a2_3 = in2(10);
a2_4 = in2(14);
a2_5 = in2(18);
a2_6 = in2(22);
a2_7 = in2(26);
a3_1 = in2(3);
a3_2 = in2(7);
a3_3 = in2(11);
a3_4 = in2(15);
a3_5 = in2(19);
a3_6 = in2(23);
a3_7 = in2(27);
a4_1 = in2(4);
a4_2 = in2(8);
a4_3 = in2(12);
a4_4 = in2(16);
a4_5 = in2(20);
a4_6 = in2(24);
a4_7 = in2(28);
dq1 = in1(6,:);
dq2 = in1(7,:);
dq3 = in1(8,:);
dq4 = in1(9,:);
dq5 = in1(10,:);
q1 = in1(1,:);
q3 = in1(3,:);
q5 = in1(5,:);
hDot = [dq1.*(a1_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(-3.0./1.6e+1)+a1_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)+a1_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)-a1_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)-a1_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a1_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a1_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a1_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a1_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)-a1_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)+a1_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)+a1_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0)+1.0)-dq5.*(a1_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a1_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a1_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a1_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a1_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a1_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a1_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a1_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a1_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a1_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a1_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a1_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))-dq3.*(a1_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)-a1_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)-a1_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)+a1_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)+a1_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)-a1_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)-a1_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)+a1_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)+a1_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)-a1_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)-a1_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1)+a1_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1));dq2-dq1.*(a2_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a2_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a2_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a2_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a2_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a2_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a2_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a2_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a2_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a2_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a2_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a2_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))-dq5.*(a2_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a2_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a2_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a2_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a2_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a2_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a2_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a2_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a2_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a2_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a2_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a2_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))-dq3.*(a2_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)-a2_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)-a2_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)+a2_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)+a2_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)-a2_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)-a2_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)+a2_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)+a2_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)-a2_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)-a2_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1)+a2_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1));-dq1.*(a3_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a3_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a3_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a3_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a3_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a3_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a3_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a3_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a3_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a3_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a3_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a3_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))-dq5.*(a3_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a3_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a3_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a3_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a3_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a3_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a3_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a3_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a3_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a3_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a3_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a3_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))+dq3.*(a3_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(-3.0./3.2e+1)+a3_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)+a3_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)-a3_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)-a3_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)+a3_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)+a3_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)-a3_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)-a3_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)+a3_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)+a3_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1)-a3_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1)+1.0);dq4-dq1.*(a4_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a4_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a4_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a4_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a4_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a4_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a4_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a4_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a4_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a4_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a4_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a4_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))-dq5.*(a4_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a4_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./1.6e+1)-a4_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a4_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./1.6e+1)+a4_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)-a4_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)-a4_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./8.0)+a4_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./8.0)+a4_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).*(1.5e+1./1.6e+1)+a4_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a4_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./1.6e+1)-a4_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./8.0))-dq3.*(a4_1.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)-a4_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^5.*(3.0./3.2e+1)-a4_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)+a4_7.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^5.*(3.0./3.2e+1)+a4_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)-a4_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)-a4_4.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^3.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^2.*(1.5e+1./1.6e+1)+a4_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^2.*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^3.*(1.5e+1./1.6e+1)+a4_5.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)-a4_6.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).*(q1.*2.0+q3+q5.*2.0+tp.*2.0).^4.*(1.5e+1./3.2e+1)-a4_2.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1)+a4_3.*1.0./(tm-tp).^6.*(q1.*2.0+q3+q5.*2.0+tm.*2.0).^4.*(q1+q3./2.0+q5+tp).*(1.5e+1./1.6e+1))];