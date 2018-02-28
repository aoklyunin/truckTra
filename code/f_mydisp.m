function f_mydisp(var,val)
disp('\begin{equation}');
disp('\tiny{');
disp(var);
try
    disp(latex(simplify(val)));
catch
     disp(val);
end
disp('}\end{equation}');
end

