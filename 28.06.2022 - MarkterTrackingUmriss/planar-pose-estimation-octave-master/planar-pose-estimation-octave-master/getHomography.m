%X, Y: 3xN Arrays: 2D points with each column representing a single point in homogeneous coordinates (x,y,w)
% retrieves original transform back given point correspondences

% (p11*X11 + p12*X21 + p13*X31)/(p31*X11 + p32*X21 + p33*X31) = Y11/Y31;
% (p21*X11 + p22*X21 + p23*X31)/(p31*X11 + p32*X21 + p33*X31) = Y21/Y31;

% p11*X11*Y31 + p12*X21*Y31 + p13*X31*Y31 + 0 + 0 + 0 - Y11*p31*X11 - Y11*p32*X21 - Y11*p33**X31 = 0 % first coloumn: x coordinate
% 0 + 0 + 0 + p21*X11*Y31 + p22*X21*Y31 + p23*X31*Y31 - Y21*p31*X11 - Y21*p32*X21 - Y21*p33**X31 = 0 % second coloumn: y coordinate

%p = [p11 ... p33]
%b = 0
%A = [X1i*Y31, X2i*Y31, X3i*Y31, 0, 0, 0, -Y1i*X1i, -Y1i*X21i, -Y1i*X3i;
%     0, 0, 0, X1i*Y31, X2i*Y31, X3i*Y31, -Y2i*X1i, -Y2i*X2i, -Y2i*X3i;
%     ...

% (see also Hartley & Zisserman Alg 3.2 page 92 in 1st edition, Alg 4.2 page 109 in 2nd edition


function H = getHomography(X, Y)

   % with i being the ith coloum (point) of X and Y
     
    numPoints = size(X, 2);
    ooo  = zeros(1,3);
    
    for i = 1:numPoints
        x = X(:,i); % access the ith point pair
        y = Y(:,i);
        A(2*i-1,:) = [x'*y(3), ooo, -x'*y(1)]; % multiplying with the w coordinate of the y vector allows w to be not 1 (which is typically the case after applying a homography)
        A(2*i,:)   = [ooo, x'*y(3), -x'*y(2)];
        
        % long form
        % A(2*i-1,:) = [X(1,i)*Y(3,i), X(2,i)*Y(3,i), X(3,i)*Y(3,i), 0, 0, 0 -Y(1,i)*X(1,i), -Y(1,i)*X(2,i), -Y(1,i)*X(3,i)];
       % A(2*i,:)   = [0, 0, 0, X(1,i)*Y(3,i), X(2,i)*Y(3,i), X(3,i)*Y(3,i), -Y(2,i)*X(1,i), -Y(2,i)*X(2,i), -Y(2,i)*X(3,i) ];
    end

    [U,D,V] = svd(A); 
    
    % check for numerical stability of result
    s = diag(D); 
    nullspace_dimension = sum(s < eps * s(1) * 1e3);
    if nullspace_dimension > 1
      fprintf('Nullspace is a bit roomy...');
    end

    % Extract homography
    H = reshape(V(:,9),3,3)'; % the solution is the eigenvector associated with the smallest eigenvalue, i.e. the last coloum of V

    % normalize so that H(3:3) is 1
    H = H(:,:) ./ H(3,3);
end
