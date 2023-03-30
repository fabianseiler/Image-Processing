clc;
clear all;


%% Image Addition
PSNR = zeros(6,1);
SSIM = zeros(6,1);
Img1 = imread('rice.png');
Img2 = imread('cameraman.tif');
Ref = imadd(Img1/2,Img2/2);
for i = (0:5)
    A = uint8(ImgAddition(uint8(Img1)/2,uint8(Img2)/2,i,8,[256 256],@AppFA));
    PSNR(i+1) = psnr(A,Ref);
    SSIM(i+1) = ssim(A,Ref);
    subplot(2,3,i+1);
    imshow(A,[]);
    title(['Appr/Exact: ', num2str(i),'/',num2str((8-i)),]);
    xlabel(['PSNR: ',num2str(psnr(A,Ref)), ', SSIM: ', num2str(ssim(A,Ref))])
end
%% Image Substraction
% Hat noch Probleme beim Scaling (vil nur halbes Img2 abziehen)
Img1 = imread('rice.png');
Img2 = imread('riceblurred.png');
Ref = imsubtract(Img1,Img2);
A = ImgSubstraction(Img1,Img2,0,8,@SSIAFA6);
imshow(Ref,[])
figure
imshow(A,[])

%% Grayscale Conversion
% Genaue Berechnung anschauen
PSNR = zeros(6,1);
SSIM = zeros(6,1);
rgb_img = imread('ngc6543a.jpg');
rgb_img = imresize(rgb_img,[256 256]);
Ref = uint8(ImgAddition(uint8(rgb_img(:,:,1))/3,ImgAddition(uint8(rgb_img(:,:,2))/3,uint8(rgb_img(:,:,3))/3,0,8,[256 256],@SSIAFA3),0,8,[256 256],@SSIAFA3));
for i = (0:5)
    A = uint8(ImgAddition(uint8(rgb_img(:,:,1))/3,ImgAddition(uint8(rgb_img(:,:,2))/3,uint8(rgb_img(:,:,3))/3,i,8,[256 256],@SSIAFA3),i,8,[256 256],@SSIAFA3));
    PSNR(i+1) = psnr(A,Ref);
    SSIM(i+1) = ssim(A,Ref);
    subplot(2,3,i+1);
    imshow(A,[]);
    title(['Appr/Exact: ', num2str(i),'/',num2str((8-i)),]);
    xlabel(['PSNR: ',num2str(psnr(A,Ref)), ', SSIM: ', num2str(ssim(A,Ref))])
end

%% Image Multiplication
Img1 = imread('cameraman.tif');
Img2 = imread('pout.tif');
Img2 = imresize(Img2, [256 256]);
A = immultiply(Img1/8,Img2/8);
imshow(A,[])

%% Quality Metrics
MED(4,8,@SSIAFA3)


%%

function Out = MED(k,n,fun)
    Out = 0;
    for i = (0:2^n-1)
        for j = (0:2^n-1)
            Out = Out + ED(i,j,k,n,fun);
        end
    end
    Out = Out/(2^(2*n));
end

function Out = ED(Int1,Int2,k,n, fun)
    Out=0;
    Ain = int2bit(Int1,n,false);
    Bin = int2bit(Int2,n,false);
    Cin = 0;
    SumFA = zeros(n,1);
    for i = [1:n]
        SumFA(i) = xor(Ain(i),xor(Bin(i),Cin));
        Cin = (Ain(i) & Bin(i)) | (Cin & (Ain(i) | Bin(i)));
    end
    SumAFA = int2bit(ApprAddition(Int1,Int2,k,n,fun),n,false);
    for j=(1:n)
        Out = Out + abs(double(SumFA(j))- double(SumAFA(j)));
    end
end

function Out = int2twoscomp(Value,n) 
% Returns -Value in 2s Complement Binary
    rest = 2^n - Value;
    Out = zeros(n,1);
    Out(1)=1;
    for i=(1:n)
        if(rest >= 2^(n-i))
            Out(i)=1;
            rest = rest - 2^(n-i);
        else
            Out(i)=0;
        end
    end    
end

function Output = ImgSubstraction(Img1,Img2,k,n, fun)
    Output = zeros(2^n,2^n);
    for i = (1:2^n)
        for j = (1:2^n)
            Output(i,j) = ApprSubstraction(Img1(i,j),Img2(i,j),k,n, fun);
        end
    end
end

function Out = ApprSubstraction(Int1,Int2,k,n, fun)
    % Calculates an n-Bit Substraction with k Approximated FA and (n-k) FA
    % k... Approximate Bits, n... Size
    Ain = int2bit(Int1,n,false);
    Bin_reverse = int2twoscomp(Int2,n);
    Bin = zeros(n,1);
    for m=(1:n)
        Bin(n-m+1) = Bin_reverse(m);
    end
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

function Output = ImgAddition(Img1,Img2,k,n,size ,fun)
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

function [Cout, Sum] = SSIAFA3(Ain,Bin,Cin)
    Cout= (Ain&Bin) | Cin;
    Sum = ~Cout;
end

function [Cout, Sum] = SSIAFA4(Ain,Bin,Cin)
    Cout= (Ain&Cin) | Bin;
    Sum = ~Cout;
end

function [Cout, Sum] = SSIAFA5(Ain,Bin,Cin)
    Cout= (Bin&Cin) | Ain;
    Sum = ~Cout;
end

function [Cout, Sum] = SSIAFA6(Ain,Bin,Cin)
    Cout= (Ain|Bin)&Cin;
    Sum = ~Cout;
end

function [Cout, Sum] = AppFA(Ain,Bin,Cin)
    % Approximate Full Adder Calculation for 1 Bit with Algorithm

    %% Alte Algorithmen
    %----Algo 1-----
    %Cout = (Ain & Bin) | (Cin & (Ain | Bin));
    %Sum = ~Cout;

    %----Algo 2----
    %Cout = (Ain|Bin) & Cin;
    %Sum = ~Cout;
    %%

    %% Neue Algorithmen
    %----Algo 3----
    %Cout= (Ain&Bin) | Cin;
    %Sum = ~Cout;

    %----Algo 4----
    %Cout= (Ain&Cin) | Bin;
    %Sum = ~Cout;

    %----Algo 5----
    %Cout= (Bin&Cin) | Ain;
    %Sum = ~Cout;

    %----Algo 6----
    %Cout= (Ain|Bin)&Cin;
    %Sum = ~Cout;
    %%

    %% Vergleiche
    %-------SIAFA1------------
    %Sum = ~Bin | (~Ain & ~Cin);
    %Cout = ~Sum;

    %-------SIAFA4------------
    Sum = ~Cin | (~Ain & ~Bin);
    Cout = ~Sum; 
    %%
end



