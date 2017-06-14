function H_out = pcafhog(H,pcasz,frame)
%将输入矩阵H经过pca进行降维
%pcasz为输出维数
    [h,w,p] = size(H);
%     pcasz = 5;%需要保留的输出维数
    H_out = zeros(h,w,pcasz);
%     if frame==58
%         disp('暂停');
%     end
%     percent=0;
    for i=1:h
        feature = reshape(H(i,:,:),w,p);
        [coed,score,latent] = pca(feature);
        coef = coed';
        H_out(i,:,:) = score*coef(:,1:pcasz);
%         pareto(100*latent/sum(latent));%调用matla画图
%         percent = percent+100*(sum(latent(1:pcasz))/sum(latent));
%         fangcha=sqrt(std(latent));
    end
%     percent = percent/h;
end