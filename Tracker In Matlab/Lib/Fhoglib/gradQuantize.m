%% gradQuantize量化O和M到O0、O1和M0、M1
function [M0,M1,O0,O1] = gradQuantize(M,O,x,h,nb,h0,norm,nOrients,full,interpolate)%量化O和M
    %gradQuantize(O+x*h,M+x*h,O0,O1,M0,M1,nb,h0,sInv2,nOrients,full,softBin>=0);
    M0=zeros(1,h);M1=zeros(1,h);O0=zeros(1,h);O1=zeros(1,h);%%
    current = (x-1)*h;
    if full
        oMult = nOrients/(2*pi);
    else
        oMult = nOrients/pi;
    end
    oMax=nOrients*nb;
    if interpolate
        for i=1:h0-3
            o=O(current+i)*oMult;o0=floor(o);od=o-o0;
            o0=o0*nb;
            if o0>=oMax,o0=0;end
            O0(i)=o0;
            o1=o0+nb;
            if o1>=oMax,o1=0;end
            O1(i)=o1;
            m=M(current+i)*norm;M1(i)=od*m;M0(i)=m-M1(i);
        end
    else
        for i=1:h0-3
            o=O(current+i)*oMult; 
			o0=floor(o+0.5); 
			o0=o0*nb;
			if o0>=oMax, o0=0;end 
			O0(i)=o0; 
			M0(i)=M(current+i)*norm; 	
			M1(i)=0.0; O1(i)=0;
        end
    end
    if interpolate
        while(1)
            if i>h0,break;end
			o=O(current+i)*oMult; 
			o0=floor(o); 
			od=o-o0;    %od幅值量化后误差
			o0=o0*nb; 
			if o0>=oMax ,o0=0;end
			O0(i)=o0;          %O = wb*hb*(取整(O*9/2PI))
			o1=o0+nb; 
			if(o1==oMax) ,o1=0;end 
			O1(i)=o1;            %O1 = wb*hb*(取整(O*9/2PI))+wb*hb
			m=M(current+i)*norm; 
			M1(i)=od*m; M0(i)=m-M1(i); %M1为梯度幅度M乘以权值后再乘以od（norm*M*od），M0=norm*M-（norm*M*od）。
            i=i+1;
        end
    else
        while(1)
            if i>h0,break;end
			o=O(current+i)*oMult;				%oMult = 9/2PI
			o0=floor(o+0.5);
			o0=o0*nb; 
			if(o0>=oMax) ,o0=0;end 
			O0(i)=o0;				%O = wb*hb*((O*9/2PI)+0.5)
			M0(i)=M(current+i)*norm; %M0为梯度幅度M乘以权值后
			M1(i)=0; O1(i)=0; %M1和O1都是0
            i=i+1;
        end
    end      
end
