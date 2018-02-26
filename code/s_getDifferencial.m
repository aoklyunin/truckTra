syms t;
syms a(t) real; syms x(t) real; syms y(t) real;
syms at real; syms xt real; syms yt real;
syms dat real; syms dxt real; syms dyt real;
syms ddat real; syms ddxt real; syms ddyt real;
syms m real; syms I real;
syms R;
syms cx real; syms cy real;

curM = [-I-m*cx^2+cy^2 -m*cy m*cx;
	m*cy -m 0;
	m*cx 0 -m;
    ];
curC =[0 dat*cx dat*cy;	
	dat*cx 0 0;
	dat*cy 0 0;
    ];

B = [zeros(3);inv(curM)];
A = [zeros(6,3) [eye(3); -inv(curM)*curC]];


disp(latex(simplify(A)));
disp(latex(simplify(B)));
