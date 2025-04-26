clear
close all

today = datestr(date(), 'yyyymmdd');
experiments = [2]; 
if any(experiments == 1) % exp 1: initialize model{{{
	steps = [1,2];
	md = runme('steps', steps);
end %}}}
if any(experiments == 2) % exp 3: inversion Weertman{{{
	steps = [3];
	costcoeffs = [1000, 1, 1e-8];
	md = runme('steps', steps, ...
		'friction', 'Weertman', ...
		'cost coefficients', costcoeffs);
end %}}}
if any(experiments == 4) % exp 4: inversion Budd{{{
	steps = [5];
	costcoeffs = [500, 1, 1e-8];
	md = runme('steps', steps, ...
		'friction', 'Budd', ...
		'cost coefficients', costcoeffs);
end %}}}

return
