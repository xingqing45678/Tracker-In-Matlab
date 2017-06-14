function [U, mu, vars] = pca_qw(X)
%PCA Principal Component Analysis
%   Supported syntaxes for tall arrays:
%
%   COEFF = pca(X) 
%   [COEFF, SCORE, LATENT] = pca(X)
%   [COEFF, SCORE, LATENT, EXPLAINED] = pca(X)
%   [COEFF, SCORE, LATENT, TSQUARED] = pca(X)
%   [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X)
%
%   Limitation:
%   Name-value pair arguments are not supported.
%
% See also PCA, PCACOV

%   Copyright 2016 The MathWorks, Inc.
d=size(X);n=d(end);d=prod(d(1:end-1));
if(~isa(X,'double')),X=double(X);end
if(n==1);mu=X;U=zeros(d,1);vars=0;return;end
mu = mean(X,ndims(X));
X = bsxfun(@minus,X,mu)/sqrt(n-1);
X = reshape(X,d,n);

m=2500;if(min(d,n)>m),X=X(:,randperm(n,m));n=m;end

if( 0 )
    [U,S]=svd(X,'econ');vars=diag(S).^2;
elseif( d>n )
    [~,SS,V]=robustSvd(X'*X);vars=diag(SS);
    U = X *V * diag(1./sqrt(vars));
else
    [~,SS,U]=robustSvd(X*X');vars=diag(SS);
end

K=vars>1e-30;vars=vars(K);U=U(:,K);
end

function [U,S,V] = robustSvd(X,trials)
    if(nargin<2),trials=100;end
    try [U,S,V] = svd(X);
    catch
        if(trials<=0),error('svd did not converge');end
        n=numel(X);j=randi(n);X(j)=X(j)+eps;
        [U,S,V]=robustSvd(X,trials-1);
    end
end
