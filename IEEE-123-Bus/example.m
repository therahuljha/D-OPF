% x = optimvar('x');
% y = optimvar('y');
% prob = optimproblem;
% prob.Objective = -x - y/3;
% prob.Constraints.cons1 = x + y <= 2;
% prob.Constraints.cons2 = x + y/4 <= 1;
% prob.Constraints.cons3 = x - y <= 2;
% prob.Constraints.cons4 = x/4 + y >= -1;
% prob.Constraints.cons5 = x + y >= 1;
% prob.Constraints.cons6 = -x + y <= 2;
% 
% sol = solve(prob)

% x = optimvar('x',2,1,'LowerBound',0);
% x3 = optimvar('x3','Type','integer','LowerBound',0,'UpperBound',1);
% prob = optimproblem;
% prob.Objective = -3*x(1) - 2*x(2) - x3;
% prob.Constraints.cons1 = x(1) + x(2) + x3 <= 7;
% prob.Constraints.cons2 = 4*x(1) + 2*x(2) + x3 == 12;
% 
% % options = optimoptions('intlinprog','Display','off');
% 
% sol = solve(prob,'intlinprog')

x = optimvar('x',3,1,'LowerBound',0);
prob = optimproblem;
prob.Objective = x(1);
cons1 = optimconstr(1);
% cons1 = 0;
for i = 1:3
    cons1 = x(i)==7;
end
%  cons1 = 7;
cons2 = 4*x(1) + 2*x(2) + x(3) == 12;
prob.Constraints.cons1 = cons1;
prob.Constraints.cons2 = cons2;

% options = optimoptions('intlinprog','Display','off');

sol = solve(prob,'intlinprog')