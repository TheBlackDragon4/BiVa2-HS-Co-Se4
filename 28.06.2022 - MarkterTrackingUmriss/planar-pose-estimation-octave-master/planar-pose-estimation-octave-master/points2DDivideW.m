function X = points2DDivideW(X)

    X(1, :) = X(1, :) ./ X(3, :);
    X(2, :) = X(2, :) ./ X(3, :);
    X(3, :) = 1;
end