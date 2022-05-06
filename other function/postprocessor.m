function [W] = postprocessor(Zn)
[uu,s,~] = svd(Zn);
s = diag(s);
r = sum(s>1e-6);
uu = uu(:, 1 : r);
s = diag(s(1 : r));
M = uu * s.^(1/2);
mm = normr(M);
rs = mm * mm';
W = rs.^(2);
end

