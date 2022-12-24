function B = mub4(x,y,z)
  % Construction of nonstandard incomplete sets of MUBs
  % The dimension is four
  % The output contains two or three bases
  % Reference: arXiv:0907.4097
  % Contact sebastien.designolle@gmail.com for questions
  B = eye(4);
  B(:,:,2) = 1/2*[1,1,1,1;1,1,-1,-1;1,-1,1i*exp(1i*x),-1i*exp(1i*x);1,-1,-1i*exp(1i*x),1i*exp(1i*x)];
  if nargin > 1
    B(:,:,3) = 1/2*[1,1,1,1;1,1,-1,-1;-exp(1i*y),exp(1i*y),exp(1i*z),-exp(1i*z);exp(1i*y),-exp(1i*y),exp(1i*z),-exp(1i*z)];
  end
end
