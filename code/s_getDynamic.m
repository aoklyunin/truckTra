
syms t;
syms a(t) real; syms x(t) real; syms y(t) real;
syms at real; syms xt real; syms yt real;
syms dat real; syms dxt real; syms dyt real;
syms ddat real; syms ddxt real; syms ddyt real;
syms m real; syms I real;
syms cx real; syms cy real;

% ������� �������� �� �� ������� � ������� �� 
A = [cos(a(t)) -sin(a(t)) x(t); 
     sin(a(t)) cos(a(t))  y(t);
     0         0          1  ];

% ������-������ ������ ���� � ������� ������� ���������
R = A*[cx;cy;1];
f_mydisp('R=',R);
% ������ ����������� ������ �������
dR = diff(R,t);

f_mydisp('dR=',dR);
f_mydisp('dR^TdR=',dR'*dR);

% ������������� �������
K = m*(dR'*dR)/2+I*diff(a(t),t)'*diff(a(t),t)/2;
f_mydisp('K=',K);

curRTR = subs(dR'*dR,[diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],[dat,dxt,dyt,at,xt,yt]);
dRTRdQ = [diff(curRTR,at);diff(curRTR,xt);diff(curRTR,yt)];
f_mydisp('dRTRdQ=',dRTRdQ);

% �������� ����������� ������������ ������� �� ���������� �����������
curK = subs(K,[diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],[dat,dxt,dyt,at,xt,yt]);
curdKdq = [diff(curK,at);diff(curK,xt);diff(curK,yt)];
dKdq = subs(curdKdq,[at,xt,yt,dat,dxt,dyt],[a(t),x(t),y(t),diff(a(t),t),diff(x(t),t),diff(y(t),t)]);
f_mydisp('dKdq=',dKdq);

% �������� ����������� ������������ ������� �� ������ ����������� ���������� ���������
curK = subs(K,[diff(a(t),t),diff(x(t),t),diff(y(t),t)],[dat,dxt,dyt]);
curdKddq = [diff(curK,dat);diff(curK,dxt);diff(curK,dyt)];
dKddq = subs(curdKddq,[dat,dxt,dyt],[diff(a(t),t),diff(x(t),t),diff(y(t),t)]);
f_mydisp('dKddq=',dKddq);


dKddqdt = diff(dKddq,t);

% �������� ����������� �������
L = dKdq - dKddqdt;

% �������� ������ ���������� ������������� ���������
v = [ddat, ddxt, ddyt, 0,0,0,0,0,0];
M = subs(L,[diff(a(t),t,t),diff(x(t),t,t),diff(y(t),t,t),diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],v);
f_mydisp('M=',M);


% �������� ������ ���������� ������������� ���������
v = [0,0,0,dat, dxt, dyt,0,0,0];
C = subs(L,[diff(a(t),t,t),diff(x(t),t,t),diff(y(t),t,t),diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],v);
f_mydisp('C=',C);

% �������� ������ ���������� ������������� ���������
v = [0,0,0,0,0,0,at,xt,yt];
G = subs(L,[diff(a(t),t,t),diff(x(t),t,t),diff(y(t),t,t),diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],v);
f_mydisp('G=',G);
