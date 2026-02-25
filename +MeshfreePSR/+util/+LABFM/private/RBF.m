function o = RBF(r, type)

switch type
    case "WendlandC2"
        o = (max(0, (1-r)).^3).*(1 + 3.*r);
    case "Gaussian"
        o = exp(-6*r.^2);
    case "Unity"
        o = ones(size(r));
    otherwise
        error("Not implemented.")
end

end