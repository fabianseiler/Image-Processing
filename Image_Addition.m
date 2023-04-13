%% Image Addition
Img1 = imread('rice.png');
Img2 = imread('cameraman.tif');

%ShowAllAlgos(Img1,Img2)0
%ShowOneAlgo(Img1,Img2,@SSIAFA1)
%[PSNR, MSSIM] = calculateMetrics(Img1,Img2,@SSIAFA4)

%% functions

function [PSNR, MSSIM] = calculateMetrics(Img1,Img2,fun)
% Calculates PSNR and MSSIM for different Approximation Degrees for one
% Algorithm
    Ref = imadd(Img1/2,Img2/2);
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,1,8,[256 256],fun));
    B = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,2,8,[256 256],fun));
    C = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,3,8,[256 256],fun));
    D = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,4,8,[256 256],fun));
    E = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,5,8,[256 256],fun));
    PSNR = [psnr(A,Ref);psnr(B,Ref);psnr(C,Ref);psnr(D,Ref);psnr(E,Ref)];
    MSSIM = [ssim(A,Ref);ssim(B,Ref);ssim(C,Ref);ssim(D,Ref);ssim(E,Ref)];
end

function ShowOneAlgo(Img1,Img2,fun)
% Shows One Algorithms with Different Approximation Degrees
    Ref = imadd(Img1/2,Img2/2);
    subplot(2,4,1);
    imshow(Img1,[]);
    title("(a)")
    subplot(2,4,2);
    imshow(Img2,[]);
    title("(b)")
    subplot(2,4,3);
    imshow(Ref,[]);
    title("(c)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,1,8,[256 256],fun));
    subplot(2,4,4);
    imshow(A,[]);
    title("(d)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,2,8,[256 256],fun));
    subplot(2,4,5);
    imshow(A,[]);
    title("(e)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,3,8,[256 256],fun));
    subplot(2,4,6);
    imshow(A,[]);
    title("(f)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,4,8,[256 256],fun));
    subplot(2,4,7);
    imshow(A,[]);
    title("(g)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,5,8,[256 256],fun));
    subplot(2,4,8);
    imshow(A,[]);
    title("(h)")
end

function ShowAllAlgos(Img1,Img2)
% Shows all the Different Image Addition Results with all Algorithms with
% an Approximation Degree of 4/8
    Ref = imadd(Img1/2,Img2/2);
    subplot(2,4,1);
    imshow(Img1,[]);
    title("(a)")
    subplot(2,4,2);
    imshow(Img2,[]);
    title("(b)")
    subplot(2,4,3);
    imshow(Ref,[]);
    title("(c)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,4,8,[256 256],@SSIAFA1));
    subplot(2,4,4);
    imshow(A,[]);
    title("(d)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,4,8,[256 256],@SSIAFA2));
    subplot(2,4,5);
    imshow(A,[]);
    title("(e)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,4,8,[256 256],@SSIAFA3));
    subplot(2,4,6);
    imshow(A,[]);
    title("(f)")
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,4,8,[256 256],@SSIAFA4));
    subplot(2,4,7);
    imshow(A,[]);
    title("(g)")
end

function void = RCA_ImageAddition()
% Image Addition with Hardcoded RCA
% Test purposes Only
PSNR = zeros(6,1);
SSIM = zeros(6,1);
Img1 = imread('rice.png');
Img2 = imread('cameraman.tif');
Ref = imadd(Img1/2,Img2/2);
A = uint8(zeros(256,256));
for i = (1:256)
    for j = (1:256)
        A(i,j) = uint8(RCA(uint8(Img1(i,j)/2),uint8(Img2(i,j)/2), @AppFA));
    end
end
imshow(A,[]);
title(['Appr/Exact: ', num2str(2),'/',num2str(8),]);
[MSSIM, ~] = ssim(A,Ref);
xlabel(['PSNR: ', num2str(psnr(A,Ref)), ', MSSIM: ' , num2str(MSSIM)])
figure
imshow(Ref,[])
end

function Output = ImgAddition(Img1,Img2,k,n,size ,fun)
% Uses ApprAddition for every Pixel in the two Grayscale Images
    Output = zeros(size(1),size(2));
    for i = (1:size(1))
        for j = (1:size(2))
            Output(i,j) = ApprAddition(Img1(i,j),Img2(i,j),k,n, fun);
        end
    end
end

function Out = ApprAddition(Int1, Int2, k, n, fun)
    % Calculates an n-Bit Addition with k Approximated FA and (n-k) FA
    % k... Approximate Bits, n... Size
    Ain = int2bit(Int1,n,false);
    Bin = int2bit(Int2,n,false);
    Cin = 0;
    Sum = zeros(n,1);
    if k>0
        for i = (1:k)
            [Cout, Sum(i)] =  fun(Ain(i),Bin(i),Cin);
            Cin = Cout;
        end
    end
    if k<8
        for j = (k+1:n)
            Sum(j) = xor(Ain(j),xor(Bin(j),Cin));
            Cin = (Ain(j) & Bin(j)) | (Cin & (Ain(j) | Bin(j)));
        end 
    end
    Out = uint8(bit2int(Sum,n,false));
end

function Out= RCA(Int1, Int2, fun)
% Testing function for RCA with stiff values
% Test purpose only
    Ain = int2bit(Int1,8,false);
    Bin = int2bit(Int2,8,false);
    Cin1 = 0;
    
    %Bit 1
    [Cin2, Sum(1)] =  fun(Ain(1),Bin(1),Cin1);

    %Bit 2
    [Cin3, Sum(2)] =  fun(Ain(2),Bin(2),Cin2);
    
    %Bit 3
    Sum(3) = xor(Ain(3),xor(Bin(3),Cin3));
    Cin4 = (Ain(3) & Bin(3)) | (Cin3 & (Ain(3) | Bin(3)));

    %Bit 4
    Sum(4) = xor(Ain(4),xor(Bin(4),Cin4));
    Cin5 = (Ain(4) & Bin(4)) | (Cin4 & (Ain(4) | Bin(4)));
     
    %Bit 5
    Sum(5) = xor(Ain(5),xor(Bin(5),Cin5));
    Cin6 = (Ain(5) & Bin(5)) | (Cin5 & (Ain(5) | Bin(5)));

    %Bit 6
    Sum(6) = xor(Ain(6),xor(Bin(6),Cin6));
    Cin7 = (Ain(6) & Bin(6)) | (Cin6 & (Ain(6) | Bin(6)));

    %Bit 7
    Sum(7) = xor(Ain(7),xor(Bin(7),Cin7));
    Cin8 = (Ain(7) & Bin(7)) | (Cin7 & (Ain(7) | Bin(7)));

    %Bit 8
    Sum(8) = xor(Ain(8),xor(Bin(8),Cin8));
    Cout = (Ain(8) & Bin(8)) | (Cin8 & (Ain(8) | Bin(8)));

    Out = uint8(bit2int(transpose(uint8(Sum)),8,false));
end

function [Cout, Sum] = SSIAFA1(Ain,Bin,Cin)
    Cout= (Ain&Bin) | Cin;
    Sum = ~Cout;
end

function [Cout, Sum] = SSIAFA2(Ain,Bin,Cin)
    Cout= (Ain&Cin) | Bin;
    Sum = ~Cout;
end

function [Cout, Sum] = SSIAFA3(Ain,Bin,Cin)
    Cout= (Bin&Cin) | Ain;
    Sum = ~Cout;
end

function [Cout, Sum] = SSIAFA4(Ain,Bin,Cin)
    Cout= (Ain|Bin)&Cin;
    Sum = ~Cout;
end

function [Cout, Sum] = AppFA(Ain,Bin,Cin)
    % Approximate Full Adder Calculation for 1 Bit with Algorithm
    %-------SIAFA1------------
    %Sum = ~Bin | (~Ain & ~Cin);
    %Cout = ~Sum;

    %-------SIAFA4------------
    Sum = ~Cin | (~Ain & ~Bin);
    Cout = ~Sum; 
    %%
end