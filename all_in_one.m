clc;
clear all;
close all;


bs=8;
x=imresize(imread('lena.jpg'),[1088,1088]);

yuv=rgb2ycbcr(x);
z=double(yuv(:,:,1));
cb=double(yuv(:,:,2));
cr=double(yuv(:,:,3));

[M, Dif, In] = Y_downsample(z);
downcb = C_downsample_A(cb);
downcbL = C_downsample_L(cb);


downcr = C_downsample_A(cr);
downcrL = C_downsample_L(cr);


sz = cat(5, M, Dif, In, downcb, downcr);
szL = cat(5, M, Dif, In, downcbL, downcrL);
temp=((bs*bs)-1);
    for i=1:temp
        [Min ,Diff, Ind, dcb, dcr] = DCT_IDCT_zig_zag(sz,i);
        [MinL ,DiffL, IndL, dcbL, dcrL] = DCT_IDCT_zig_zag(szL,i);
        
        
        upz = Y_interpolation(Min, Diff, Ind);
        
        upcb = C_Bi_linear_interpolation(dcb);
        upcbL = C_Bi_linear_interpolation(dcbL);
        %upcbR = C_Bi_linear_interpolation(dcbR);
        upcbt = upsampling(dcb);
        upcbt1 = upsampling(dcbL);
        
        upcr = C_Bi_linear_interpolation(dcr);
        upcrL = C_Bi_linear_interpolation(dcrL);
        %upcrR = C_Bi_linear_interpolation(dcrR);
        upcrt = upsampling(dcr);
        upcrt2 = upsampling(dcrL);
                
        RGB = yuv2rgb(upz,upcb,upcr);
        RGBL = yuv2rgb(upz,upcbL,upcrL);
        RGBR = yuv2rgb(z,upcbt,upcrt);
        RGBR1 = yuv2rgb(z,upcbt1,upcrt2);
        [psnr_d, ms_d]=measerr(x,RGB);
        [psnr_dL, ms_dL]=measerr(x,RGBL);
        [psnr_dR, ms_dR]=measerr(x,RGBR);
        [psnr_dR1, ms_dR1]=measerr(x,RGBR1);
        
        PSNR(i) =   psnr_d;
        PSNRL(i) =   psnr_dL;
        PSNRRi(i) =   psnr_dR;
        PSNRR1(i) =   psnr_dR1;
        
        
    end
    
    
result_analysis(PSNR,PSNRL,PSNRRi,PSNRR1);%,MSE);  

 
 
  





function img = C_downsample_A(q)
[r c]=size(q);
ra = r-1;
ca = c-1;

i=1;
for l=1:2:ra
    j=1;
    for k=1:2:ca
        img(i,j)=round((q(l,k)+q(l,k+1)+q(l+1,k)+q(l+1,k+1))/4,0);
        j=j+1;
    end
    i=i+1;
    
end
end

function img = C_downsample_L(q)
[r c]=size(q);
ra = r-1;
ca = c-1;

i=1;
for l=1:2:ra
    j=1;
    for k=1:2:ca
        img(i,j)=round((q(l,k)+q(l+1,k))/2,0);
        j=j+1;
    end
    i=i+1;
    
end
end

function img = C_downsample_R(q)
[r c]=size(q);
ra = r-1;
ca = c-1;

i=1;
for l=1:2:ra
    j=1;
    for k=1:2:ca
        img(i,j)=round((q(l,k+1)+q(l+1,k))/2,0);
        j=j+1;
    end
    i=i+1;
    
end
end


function img2 = upsampling(w)
[r c]=size(w);
ra = r*2;
ca = c*2;
i=1;

for l=1:2:ra
    j=1;
    for k=1:2:ca
        temp=w(i,j);
        img2(l,k)=temp;
        img2(l,k+1)=temp;
        img2(l+1,k)=temp;
        img2(l+1,k+1)=temp;
        
        j=j+1;
    end
    i=i+1;
    
end
end

function img1 = C_Bi_linear_interpolation(w)
[ro co]=size(w);
%r = ro ;
%c = co ;
ra = ro*2;
ca = co*2;
i=1;

for l=1:2:ra
    j=1;
    for k=1:2:ca
        if i == ro | i == 1 | j == co | j == 1
            temp=w(i,j);
            img1(l,k)=temp;
            img1(l,k+1)=temp;
            img1(l+1,k)=temp;
            img1(l+1,k+1)=temp;
        else
            r = i;
            c = j;
            
            u1 = 0.1875*(w(r, c-1)) + 0.0625*(w(r-1, c-1)) + 0.1875*(w(r-1, c));
            u2 = 0.1875*(w(r, c+1)) + 0.0625*(w(r-1, c+1)) + 0.1875*(w(r-1, c));
            u3 = 0.1875*(w(r, c-1)) + 0.0625*(w(r-1, c-1)) + 0.1875*(w(r+1, c));
            u4 = 0.1875*(w(r, c+1)) + 0.0625*(w(r+1, c+1)) + 0.1875*(w(r+1, c));
            
            img1(l,k)= 0.5625*(w(r, c)) + u1;
            img1(l,k+1)= 0.5625*(w(r, c)) + u2;
            img1(l+1,k)= 0.5625*(w(r, c)) + u3;
            img1(l+1,k+1)= 0.5625*(w(r, c)) + u4;
            
            
            
        end
        
        j=j+1;
    end
    i=i+1;
    
end
end

function [Max, diff, q] = Y_downsample(z)
[r c]=size(z);
ra = r-1;
ca = c-1;

i=1;
for l=1:2:ra
    j=1;
    for k=1:2:ca
        t = [z(l,k) z(l,k+1) z(l+1,k) z(l+1,k+1)];
        [Min(i,j),Mini(i,j)] = min(t);
        [Max(i,j),Maxi(i,j)] = max(t);
        [s, in] = sort(t);
        if in == [1 2 3 4]
            q(i,j) = 1;
        elseif in == [1 2 4 3]
            q(i,j) = 2;
        elseif in == [1 3 2 4]
            q(i,j) = 3;
        elseif in == [1 3 4 2]
            q(i,j) = 4;
        elseif in == [1 4 2 3]
            q(i,j) = 5;
        elseif in == [1 4 3 2]
            q(i,j) = 6;
        elseif in == [2 1 3 4]
            q(i,j) = 7;
        elseif in == [2 1 4 3]
            q(i,j) = 8;
        elseif in == [2 3 1 4]
            q(i,j) = 9;
        elseif in == [2 3 4 1]
            q(i,j) = 10;
        elseif in == [2 4 1 3]
            q(i,j) = 11;
        elseif in == [2 4 3 1]
            q(i,j) = 12;
        elseif in == [3 1 2 4]
            q(i,j) = 13;
        elseif in == [3 1 4 2]
            q(i,j) = 14;
        elseif in == [3 2 1 4]
            q(i,j) = 15;
        elseif in == [3 2 4 1]
            q(i,j) = 16;
        elseif in == [3 4 1 2]
            q(i,j) = 17;
        elseif in == [3 4 2 1]
            q(i,j) = 18;
        elseif in == [4 1 2 3]
            q(i,j) = 19;
        elseif in == [4 1 3 2]
            q(i,j) = 20;
        elseif in == [4 2 1 3]
            q(i,j) = 21;
        elseif in == [4 2 3 1]
            q(i,j) = 22;
        elseif in == [4 3 1 2]
            q(i,j) = 23;
        elseif in == [4 3 2 1]
            q(i,j) = 24;
        end
        j=j+1;
    end
    i=i+1;
    
end

diff = Max - Min ;
end

function im = Y_interpolation(max, dif, ind)

v = round(dif/4);
M = max - dif;
nxt1 = M + v;
nxt2 = max - v;

[r c]=size(dif);
ra1 = r*2;
ca1 = c*2;
ra = ra1;
ca = ca1;
i=1;
for l=1:2:ra
    j = 1;
    for k=1:2:ca 
        if ind(i,j) ==1
            a = [1 2 3 4];
        elseif ind(i,j) ==2
            a = [1 2 4 3];
        elseif ind(i,j) ==3
            a = [1 3 2 4];
        elseif ind(i,j) ==4
            a = [1 3 4 2];
        elseif ind(i,j) ==5
            a = [1 4 2 3];
        elseif ind(i,j) ==6
            a = [1 4 3 2];
        elseif ind(i,j) ==7
            a = [2 1 3 4];
        elseif ind(i,j) ==8
            a = [2 1 4 3];
        elseif ind(i,j) ==9
            a = [2 3 1 4];
        elseif ind(i,j) ==10
            a = [2 3 4 1];
        elseif ind(i,j) ==11
            a = [2 4 1 3];
        elseif ind(i,j) ==12
            a = [2 4 3 1];
        elseif ind(i,j) ==13
            a = [3 1 2 4];
        elseif ind(i,j) ==14
            a = [3 1 4 2];
        elseif ind(i,j) ==15
            a = [3 2 1 4];
        elseif ind(i,j) ==16
            a = [3 2 4 1];
        elseif ind(i,j) ==17
            a = [3 4 1 2];
        elseif ind(i,j) ==18
            a = [3 4 2 1];
        elseif ind(i,j) ==19
            a = [4 1 2 3];
        elseif ind(i,j) ==20
            a = [4 1 3 2];
        elseif ind(i,j) ==21
            a = [4 2 1 3];
        elseif ind(i,j) ==22
            a = [4 2 3 1];
        elseif ind(i,j) ==23
            a = [4 3 1 2];
        elseif ind(i,j) ==24
            a = [4 3 2 1];
        elseif ind(i,j) <=0
            a = [1 2 3 4];
        elseif ind(i,j) >=25
            a = [4 3 2 1];
        else
            a = [4 3 2 1];
        end
        for n=1:4
            if n==1
                if a(n)==1
                    im(l,k)=M(i,j);
                elseif a(n)==2
                    im(l,k+1)=M(i,j);
                elseif a(n)==3
                    im(l+1,k)=M(i,j);
                elseif a(n)==4
                    im(l+1,k+1)=M(i,j);
                end
            elseif n==2
                if a(n)==1
                    im(l,k)=nxt1(i,j);
                elseif a(n)==2
                    im(l,k+1)=nxt1(i,j);
                elseif a(n)==3
                    im(l+1,k)=nxt1(i,j);
                elseif a(n)==4
                    im(l+1,k+1)=nxt1(i,j);
                end
            elseif n==3
                if a(n)==1
                    im(l,k)=nxt2(i,j);
                elseif a(n)==2
                    im(l,k+1)=nxt2(i,j);
                elseif a(n)==3
                    im(l+1,k)=nxt2(i,j);
                elseif a(n)==4
                    im(l+1,k+1)=nxt2(i,j);
                end
            elseif n==4
                if a(n)==1
                    im(l,k)=max(i,j);
                elseif a(n)==2
                    im(l,k+1)=max(i,j);
                elseif a(n)==3
                    im(l+1,k)=max(i,j);
                elseif a(n)==4
                    im(l+1,k+1)=max(i,j);
                end        
            end
        end
        j = j+1;
    end
    i = i+1;
end

end


function RGB = yuv2rgb(upz,upcb,upcr)

    %R = 1.164*(y -16)+ 1.596* (cr -128);
    %G = 1.164*(y -16)- 0.391* (cb -128)- 0.813*(cr -128);
    %B = 1.164*(y -16)+ 2.018* (cb -128);
    f = cat(3, upz,upcb,upcr);
    g = uint8(f);
    RGB = ycbcr2rgb(g);
    %g = cat(3, R, G, B);
    %RGB = uint8(g);

end

function [Mx ,Di, Ix, dcb, dcr] = DCT_IDCT_zig_zag(sz,i)

    
    bs=8;
    temp=((bs*bs)-1);
    g = double(sz(:,:,1));
    g1 = double(sz(:,:,2));
    g2 = double(sz(:,:,3));
    g3 = double(sz(:,:,4));
    g4 = double(sz(:,:,5));
    [ra, ca]=size(g);
    cr1=mat2cell(g,(bs*ones(1,ra/bs)),(bs*ones(1,ca/bs)));
    cr2=mat2cell(g1,[bs*ones(1,ra/bs)],[bs*ones(1,ca/bs)]);
    cr3=mat2cell(g2,[bs*ones(1,ra/bs)],[bs*ones(1,ca/bs)]);
    cr4=mat2cell(g3,[bs*ones(1,ra/bs)],[bs*ones(1,ca/bs)]);
    cr5=mat2cell(g4,[bs*ones(1,ra/bs)],[bs*ones(1,ca/bs)]);
    c=[1 2 3 2 1 1 2 3 4 5 4 3 2 1 1 2 3 4 5 6 7 6 5 4 3 2 1 1 2 3 4 5 6 7 8 ... 
        8 7 6 5 4 3 2 3 4 5 6 7 8 8 7 6 5 4 5 6 7 8 8 7 6 7 8 8 ];
    d=[2 1 1 2 3 4 3 2 1 1 2 3 4 5 6 5 4 3 2 1 1 2 3 4 5 6 7 8 7 6 5 4 3 2 1 ....
        2 3 4 5 6 7 8 8 7 6 5 4 3 4 5 6 7 8 8 7 6 5 6 7 8 8 7 8 ];
    
    
        for m=1:ra/bs
            for n=1:ca/bs
                data_in1=cr1{m,n};
                data_in2=cr2{m,n};
                data_in3=cr3{m,n};
                data_in4=cr4{m,n};
                data_in5=cr5{m,n};
                dct_fwd1=dct2(data_in1);
                dct_fwd2=dct2(data_in2);
                dct_fwd3=dct2(data_in3);
                dct_fwd4=dct2(data_in4);
                dct_fwd5=dct2(data_in5);
                    for j=i:1:temp
                        dct_fwd1(c(j),d(j))=0;
                        dct_fwd2(c(j),d(j))=0;
                        for jj=40:1:temp
                            dct_fwd3(c(j),d(j))=0;
                        end
                        dct_fwd4(c(j),d(j))=0;
                        dct_fwd5(c(j),d(j))=0;
                    end
        		dct_inv1=idct2((dct_fwd1));
                dct_inv2=idct2((dct_fwd2));
                dct_inv3=idct2((dct_fwd3));
                dct_inv4=idct2((dct_fwd4));
                dct_inv5=idct2((dct_fwd5));
        		data_out_d1{m,n}=dct_inv1;
                data_out_d2{m,n}=dct_inv2;
                data_out_d3{m,n}=dct_inv3;
                data_out_d4{m,n}=dct_inv4;
                data_out_d5{m,n}=dct_inv5;
            end
        end
        Mx=uint8(cell2mat(data_out_d1));
        Di=uint8(cell2mat(data_out_d2));
        Ix=uint8(cell2mat(data_out_d3));
        dcb=uint8(cell2mat(data_out_d4));
        dcr=uint8(cell2mat(data_out_d5));
        
end

function result_analysis(PSNR,PSNRL,PSNRR,PSNRR1)
    
    c = length(PSNR);
    f = 0:c-1;
    figure
    
    plot(f,PSNR,'--r',f,PSNRL,'--b',f,PSNRR,'--g',f,PSNRR1,'--k');
    xlabel('number of values');
    ylabel('PSNR value');
    title('PSNR for various methods');
    legend('PSNR ALL with luma modification','PSNR LEFT with luma modification','PSNR All without luma modification','PANR LEFT without luma modification');
    %plot(PSNRL,'--r','Linewidth',2),xlabel('Number of transform coefficents retained'),ylabel('PSNR');
    
   
    %plot(MSE,'--r','Linewidth',2),xlabel('Number of transform coefficents retained'),ylabel('MSE');
    
end


