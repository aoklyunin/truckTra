syms t;
syms a(t) real; syms x(t) real; syms y(t) real;
syms at real; syms xt real; syms yt real;
syms dat real; syms dxt real; syms dyt real;
syms ddat real; syms ddxt real; syms ddyt real;
syms m real; syms I real;
syms R;
syms cx real; syms cy real;

curM = 	[m*(cx^2+cy^2)-I  m*(cx*sin(at)-cy*cos(at)) m*(cx*cos(at)-cy*sin(at));
	m*(-cx*sin(at)-cy*cos(at)) m 0;
	m*(cx*cos(at)-cy*sin(at)) 0 m];
    
curC =[0 (dat+1)*cx*cos(at)+(dat-1)*cy*sin(at) (dat+1)*(cx*sin(at)+cy*cos(at));
	dat*(cy*sin(at)-cx*cos(at)) 0 0;
	-dat*(cx*sin(at)+cy*cos(at)) 0 0
    ];

f_mydisp('M=',curM);
f_mydisp('C=',curC);

f_mydisp('inv(curM)=',inv(curM));
D = inv(curM)*curC;
f_mydisp('inv(curM)*curC=',inv(curM)*curC);
simplifY(curC)
