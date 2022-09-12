% This MATLAB-function was used to obtain the results of the article "Arithmetic Monodromy
% in Sp(2n)" by Jitendra Bajpai, Daniele Dona and Martin Nitsche.
% It searches in a given matrix group Gamma, Zariski-dense in Sp(2n), for an element gamma,
% such that T=1+vR*vL and its gamma-conjugate satisfy the preconditions of Lemma 1.
%
% Arguments:
%   - A, B   are matrices generating the group Gamm
%   - vL, vR are row/column matrices such that T:=1+vR*vL is an element of Gamma
%   - bnd    is an integer upper bound (e.g. 1000).
%            The search is pruned wherever gamma*vR has entries of absolute values >= bnd.
%   - w      is an element of Gamma, the starting point for the search (e.g. identity-matrix)
% Return values:
%   - found  is a  boolean indicating success
%   - its    is an integer denoting the depth at which the search stopped
%   The element gamma itself is not returned, but it can be tracked down by running the
%   function repeatedly, iteratively decreasing the its-value, by supplying varying w
%
% (Note: The run time could be further improved by using Bruno Luong's "Merge sorted arrays")

function [found, its] = witnessSearch(A, B, vL, vR, bnd, w)
  [allR newR] = deal(zeros(0,size(A,1)), w*vR);
  [found its] = deal(false, 0);
  while size(newR, 2) > 0
    if rank([vR newR(:,abs(vL*newR) == 0)]) > 1
      found = true;
      return
    end
    newR = [A*newR B*newR inv(A)*newR inv(B)*newR];
    newR = unique(newR(:, all(abs(newR) < bnd, 1))', 'rows');
    [allR idx] = sortrows([newR; allR]);
    mask = [any(allR(1:size(allR, 1)-1, :) ~= allR(2:size(allR, 1), :), 2);true];
    [allR newR its] = deal(allR(mask, :), allR(mask & idx<=size(newR,1), :)', its+1);
    display([its, size(newR, 2), size(allR, 1)]);
  end
end 
