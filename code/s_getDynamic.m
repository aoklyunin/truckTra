
syms t;
syms a(t) real; syms x(t) real; syms y(t) real;
syms at real; syms xt real; syms yt real;
syms dat real; syms dxt real; syms dyt real;
syms ddat real; syms ddxt real; syms ddyt real;
syms m real; syms I real;
syms cx real; syms cy real;

% матрица перехода из СК тележки в базовую СК 
A = [cos(a(t)) -sin(a(t)) x(t); 
     sin(a(t)) cos(a(t))  y(t);
     0         0          1  ];

% радиус-вектор центра масс в базовой системе координат
R = A*[cx;cy;1];
% первая производная радиус вектора
dR = diff(R,t);

% потенциальная энергия
K = m*(dR'*dR)/2+I*diff(a(t),t)'*diff(a(t),t)/2;
disp(simplify(K));

% получаем производную кинетической энергии по обобщённым координатам
curK = subs(K,[a(t),x(t),y(t)],[at,xt,yt]);
curdKdq = [diff(K,at);diff(K,xt);diff(K,yt)];
dKdq = subs(curdKdq,[at,xt,yt],[a(t),x(t),y(t)]);
disp(latex(simplify(dKdq)));

% получаем производную кинетической энергии по первой производной обобщённых координат
curK = subs(K,[diff(a(t),t),diff(x(t),t),diff(y(t),t)],[dat,dxt,dyt]);
curdKddq = [diff(curK,dat);diff(curK,dxt);diff(curK,dyt)];
dKddq = subs(curdKddq,[dat,dxt,dyt],[diff(a(t),t),diff(x(t),t),diff(y(t),t)]);
disp(latex(simplify(dKddq)));

dKddqdt = diff(dKddq,t);

% получаем Лагшранжиан системы
L = dKdq - dKddqdt;

% получаем первую компоненту динамического уравнения
v = [ddat, ddxt, ddyt, 0,0,0,0,0,0];
M = subs(L,[diff(a(t),t,t),diff(x(t),t,t),diff(y(t),t,t),diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],v);
disp(latex(simplify(M)));

% получаем третью компоненту динамического уравнения
v = [0,0,0,dat, dxt, dyt,0,0,0];
C = subs(L,[diff(a(t),t,t),diff(x(t),t,t),diff(y(t),t,t),diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],v);
disp(latex(simplify(C)));

% получаем первую компоненту динамического уравнения
v = [0,0,0,0,0,0,at,xt,yt];
G = subs(L,[diff(a(t),t,t),diff(x(t),t,t),diff(y(t),t,t),diff(a(t),t),diff(x(t),t),diff(y(t),t),a(t),x(t),y(t)],v);
disp(latex(simplify(G)));

