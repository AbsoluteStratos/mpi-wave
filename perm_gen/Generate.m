% Gaussian random fields with exponential correlation function
% Args: l = length scale, nsl = KLE discrization, nimg = num pixles in
% image
function y1 = Generate(l, nsl, nimg)
%==========================================================================
randn('seed', 1e5);
% Get input info:
x_surf = repmat(linspace(0, 1, nsl)', 1, nsl);
y_surf = x_surf';
x = [x_surf(:), y_surf(:)];
ns = size(x, 1);
r = sqrt(dist2(x, x));
% Get basis functions
c = exp(-r / l);
% eigs returns results in ASCENDING order
[e_vec, e_val] = eig(c);
e_val = diag(e_val);
e_val = e_val(end : -1 : 1)';
e_vec_trans = e_vec(:, end : -1 : 1)';
n_basis_functions = ns;
e_val_truncated = e_val(1 : n_basis_functions);
% NB: ROWS of basis are the functions
basis = e_vec_trans(1 : n_basis_functions, :);
coeffs = sqrt(e_val_truncated);
% Batching:
batch_bounds = [1,1];
ns = nsl^2;
n = 1;
transform = @(x) x.^2 + 0.2;
omega = randn(n, ns);
y_vec = (coeffs .* omega) * basis;

% Rescale to 0-0.1
y0 = 0.1*rescale(y_vec);
ySmall = reshape(y0', [nsl, nsl]);
y = imresize(ySmall,(nimg/nsl),'bilinear');
contourf(y);
colorbar;
y1 = transform(y);
size(y)
fileID = fopen('mat.bin','w');
fwrite(fileID, y, 'double');
fclose(fileID);

end



