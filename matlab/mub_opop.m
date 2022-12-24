function B = mub_opop(d)
  % Construction of nonstandard complete sets of MUBs
  % The dimension d is an odd power of an odd prime
  % The output contains d+1 bases
  % Reference: example 3.5b of arXiv:1104.3370
  % Contact sebastien.designolle@gmail.com for questions
  f = factor(d);
  p = f(1);
  if ~mod(p,2)
    error('p should be odd')
  end
  r = find(f-p*ones(size(f)),1)-1;
  if ~isempty(r)
    error('d should not be composite')
  end
  if ~mod(r,2)
    error('r should be odd')
  end
  r = length(f);
  B = zeros(d,d,d+1,(r+1)/2);
  for s = 0:(r-1)/2
    B(:,:,1,s+1) = eye(d);
  end
  gamma = exp(2*1i*pi/p);
  % http://fr.mathworks.com/matlabcentral/fileexchange/32872-a-toolbox-for-simple-finite-field-operation
  g = gf(p,r);
  for x = 0:d-1
    for a = 0:d-1
      for l = 0:d-1
        for s = 0:(r-1)/2
          B(:,a+1,x+2,s+1) = B(:,a+1,x+2,s+1) ...
            +gamma^(double(g.tr(g.mult(a,l))))...
            *gamma^(double(g.tr(g.add(g.mult(x,g.power(l,p^(r-s)+1)),g.mult(g.power(x,p^s),g.power(l,p^s+1)))))*(p+1)/2)*B(:,l+1,1)/sqrt(d);
        end
      end
    end
  end
end
