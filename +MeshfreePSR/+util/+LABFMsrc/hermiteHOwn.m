function f = hermiteHOwn(order, x)

switch order
    case 0
        f = ones(size(x));
    case 1
        f = 2*x;
    case 2
        f = 4*x.^2 - 2;
    case 3
        f = 8*x.^3 - 12*x;
    otherwise
        error("Not implemtend");
end


end