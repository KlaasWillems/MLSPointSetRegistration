function Xt = getMonomials(dxjiVec, dyjiVec, order)

switch order
    case 1
        Xt = [dxjiVec'; dyjiVec'];
    case 2
        Xt = [dxjiVec'; dyjiVec'; dxjiVec' .* dxjiVec' /2; dxjiVec' .* dyjiVec'; dyjiVec' .* dyjiVec' /2];
    case 3
        Xt = [  dxjiVec'; dyjiVec'; 
                dxjiVec' .* dxjiVec' /2;
                dxjiVec' .* dyjiVec'; 
                dyjiVec' .* dyjiVec' /2;
                dxjiVec' .* dxjiVec' .* dxjiVec' / 6;
                dxjiVec' .* dxjiVec' .* dyjiVec' / 2;
                dxjiVec' .* dyjiVec' .* dyjiVec' / 2;
                dyjiVec' .* dyjiVec' .* dyjiVec' / 6
                ];
    otherwise 
        error("Not implemented.")
end

end