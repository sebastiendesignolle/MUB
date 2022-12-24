function B = mub(d)
  % Construction of the standard complete set of MUBs
  % The dimension d can be any integer greater than two
  % The output contains min_i p_i^r_i+1 bases where d = p_1^r_1*...*p_n^r_n
  % Reference: arXiv:1004.3348
  % Contact sebastien.designolle@gmail.com for questions
  if d < 2
    error('p should be greater than two')
  end
  f = factor(d);
  p = f(1);
  r = find(f-p*ones(size(f)),1)-1;
  if ~isempty(r)
    B_aux1 = mub(p^r);
    B_aux2 = mub(d/p^r);
    k = min(size(B_aux1,3),size(B_aux2,3));
    B = zeros(d,d,k);
    for j = 1:k
      B(:,:,j) = kron(B_aux1(:,:,j),B_aux2(:,:,j));
    end
  else
    r = length(f);
    B = zeros(d,d,d+1);
    B(:,:,1) = eye(d);
    gamma = exp(2*1i*pi/p);
    % http://fr.mathworks.com/matlabcentral/fileexchange/32872-a-toolbox-for-simple-finite-field-operation
    g = gf(p,r);
    if p == 2
      for i = 1:d
        for k = 0:d-1
          for q = 0:d-1
            aux = 1;
            q_bin = flip(dec2bin(q,r)-'0');
            for m = 0:r-1
              for n = 0:r-1
                aux = aux*conj(1i^(double(g.mult(i-1,g.mult(q_bin(m+1)*2^m,q_bin(n+1)*2^n)))));
              end
            end
            B(:,k+1,i+1) = B(:,k+1,i+1)+gamma^(double(g.mult(q,k)))*aux*B(:,q+1,1)/sqrt(d);
          end
        end
      end
    else
      for i = 1:d
        for k = 0:d-1
          for q = 0:d-1
            B(:,k+1,i+1) = B(:,k+1,i+1)+gamma^(double(g.sub(0,g.mult(q,k))))*gamma^(double(g.mult(g.mult(i-1,g.mult(q,q)),g.inv(2))))*B(:,q+1,1)/sqrt(d);
          end
        end
      end
    end
  end
