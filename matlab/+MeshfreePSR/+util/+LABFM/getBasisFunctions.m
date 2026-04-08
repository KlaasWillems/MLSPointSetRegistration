function Wt = getBasisFunctions(dxjiVec, dyjiVec, hi, order)
% Returns the vector Wji containing the LABFM basis functions

sqrt2 = sqrt(2);
xEval = dxjiVec/(hi*sqrt2);
yEval = dyjiVec/(hi*sqrt2);

switch order
    case 1
        % Xji = [xji; yji];
        v1 = hermiteHOwn(1, xEval).*hermiteHOwn(0, yEval)/sqrt2;
        v2 = hermiteHOwn(0, xEval).*hermiteHOwn(1, yEval)/sqrt2;
        Wt = [v1, v2];
    case 2
        % Xji = [xji; yji; xji*xji/2; xji*yji; yji*yji/2];
        v1 = hermiteHOwn(1, xEval).*hermiteHOwn(0, yEval)/sqrt2;
        v2 = hermiteHOwn(0, xEval).*hermiteHOwn(1, yEval)/sqrt2;
        v3 = hermiteHOwn(2, xEval).*hermiteHOwn(0, yEval)/2;
        v4 = hermiteHOwn(1, xEval).*hermiteHOwn(1, yEval)/2;
        v5 = hermiteHOwn(0, xEval).*hermiteHOwn(2, yEval)/2;
        Wt = [v1, v2, v3, v4, v5];
    case 3
        % Xji = [xji; yji; xji*xji/2; xji*yji; yji*yji/2; xji*xji*xji/6; xji*xji*yji/2; xji*yji*yji/2; yji*yji*yji/6];
        denum = sqrt(2^3);
        v1 = hermiteHOwn(1, xEval).*hermiteHOwn(0, yEval)/sqrt2;
        v2 = hermiteHOwn(0, xEval).*hermiteHOwn(1, yEval)/sqrt2;
        v3 = hermiteHOwn(2, xEval).*hermiteHOwn(0, yEval)/2;
        v4 = hermiteHOwn(1, xEval).*hermiteHOwn(1, yEval)/2;
        v5 = hermiteHOwn(0, xEval).*hermiteHOwn(2, yEval)/2;
        v6 = hermiteHOwn(3, xEval).*hermiteHOwn(0, yEval)/denum;
        v7 = hermiteHOwn(2, xEval).*hermiteHOwn(1, yEval)/denum;
        v8 = hermiteHOwn(1, xEval).*hermiteHOwn(2, yEval)/denum;
        v9 = hermiteHOwn(0, xEval).*hermiteHOwn(3, yEval)/denum;
        Wt = [v1, v2, v3, v4, v5, v6, v7, v8, v9];
    otherwise 
        error("Not implemented.")
end

rEval = sqrt(dxjiVec .* dxjiVec + dyjiVec .* dyjiVec)/(hi);
rbfji = RBF(rEval, "Gaussian");
Wt = Wt .* rbfji;

end