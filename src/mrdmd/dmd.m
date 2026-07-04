function [modes, D, b] = dmd(X, r)
%DMD Dynamic mode decomposition with fixed truncation rank.

if nargin < 2 || isempty(r)
    r = 25;
end

X1 = X(:, 1:end-1);
Y = X(:, 2:end);
r = min(r, min(size(X1)));
[U, S, V] = svds(X1, r);
%{
thresh = 1e-8 * sing_vals(1);
r = sum(sing_vals > thresh);
r = min([25, r, size(U, 2)]);
%}
if r == 0, r = 1; end
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);

A = (U' * Y * V) / S;
[W, D] = eig(A);

modes = single(Y * V * (S \ W));
b = single(pinv(modes) * X1(:, 1));

end
