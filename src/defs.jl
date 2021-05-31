# Defs

export linterp
function linterp(x1, y1, x2, y2, x)
    R = typeof(zero(y1) + zero(y2))
    return R(x2 - x) ./ R(x2 - x1) .* y1 + R(x - x1) ./ R(x2 - x1) .* y2
end
