function res = is_mu(B)
  % Check whether the input is indeed mutually unbiased
  tol = 1e-10;
  d = size(B,1);
  k = size(B,3);
  res = true;
  for i = 1:k
    res = res && all(all(abs(B(:,:,i)'*B(:,:,i) - eye(d)) < tol));
    for j = (i+1):k
      res = res && all(all(abs(sqrt(d)*abs(B(:,:,i)'*B(:,:,j)) - ones(d)) < tol));
    end
  end
end
