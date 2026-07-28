function [c, ceq] = LyapunovConstraints(x, F_partial)
%%
% This function is used to define the Lyapunov constraints.
%%

c = F_partial(x) + 1*[0.001, 0.001, 0.001];    % c(x) <= 0               
ceq = [];                               % for Equality constraint, and [] for no Equality constraint.

%%
