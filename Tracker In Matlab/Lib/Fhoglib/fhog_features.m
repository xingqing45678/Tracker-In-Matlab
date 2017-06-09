%% compute FHOG features
function H = fhog_features(I,M,O,binSize,nOrients,softBin,useHog,clip)
% INPUTS
%  M        - [hxw] gradient magnitude at each location (see gradientMag.m)
%  O        - [hxw] gradient orientation in range defined by param flag
%  binSize  - [8] spatial bin size
%  nOrients - [9] number of orientation bins
%  softBin  - [1] set soft binning (odd: spatial=soft, >=0: orient=soft)
%  useHog   - [0] 1: compute HOG (see hog.m), 2: compute FHOG (see fhog.m)
%  clipHog  - [.2] value at which to clip hog histogram bins
%  full     - [false] if true expects angles in [0,2*pi) else in [0,pi)
%
% OUTPUTS
%  H        - [w/binSize x h/binSize x nOrients] gradient histograms
    if(nargin<3 || isempty(I)||isempty(M)||isempty(O)),  return; end
    if(nargin<4 || isempty(binSize)), binSize=8; end
    if(nargin<5 || isempty(nOrients)), nOrients=9; end
    if(nargin<6 || isempty(softBin)), softBin=-1; end
    if(nargin<7 || isempty(useHog)),useHog=2;end
    if(nargin<8 || isempty(clip)),clip=0.2;end
    %compute unnormalized constrast sensitive histograms	计算非归一化对比度敏感的直方图
    [h,w] = size(I);hb = floor(h/binSize);wb = floor(w/binSize);nb = hb*wb;nbo = nb*nOrients;
    R2 = zeros(1,wb*hb*nOrients);H1=zeros(1,wb*hb*(nOrients * 3 + 5));H=zeros(hb,wb,(nOrients * 3 + 5));
    R1 = gradientHist(M,O,h,w,binSize,nOrients*2,softBin,1);%调用梯度直方图函数
    %compute unnormalized contrast insensitive histograms	计算非归一化对比度不敏感的直方图
    for o=1:nOrients
        for x=1:nb
            R2((o-1)*nb+x)=R1((o-1)*nb+x)+R1(((o-1)+nOrients)*nb+x);
        end
    end
    %compute block normalization values 计算块归一化值
    N = hogNormMatrix(R2, nOrients, hb, wb, binSize);
    %normalized histograms and texture channels 归一化直方图和纹理通道
	H1=hogChannels(H1,nbo,0, R1, N, hb, wb, nOrients * 2, clip, 1);
	H1=hogChannels(H1,nbo,2, R2, N, hb, wb, nOrients * 1, clip, 1);
	H1=hogChannels(H1,nbo,3, R1, N, hb, wb, nOrients * 2, clip, 2);
    for i=1:nOrients * 3 + 5
        for j=1:wb
            for k=1:hb
                H(k,j,i)=H1((i-1)*hb*wb+(j-1)*hb+k);
            end
        end
    end
end



