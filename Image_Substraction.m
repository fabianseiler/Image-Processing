%% Image Substraction
Img1 = imread('images\motion02.512.tiff');
Img2 = imread('images\motion05.512.tiff');

%ShowAllAlgos(Img1,Img2)
%ShowOneAlgo(Img1,Img2,@SSIAFA4)
%[PSNR, MSSIM] = calculateMetrics(Img1,Img2,@SSIAFA4)

%% functions
function ShowAllAlgos(Img1,Img2)
% Shows all the Different Image Addition Results with all Algorithms with
% an Approximation Degree of 4/8
    Ref = imsubtract(Img1,Img2);
    subplot(2,4,1);
    imshow(Img1,[]);
    title("(a)")
    subplot(2,4,2);
    imshow(Img2,[]);
    title("(b)")
    subplot(2,4,3);
    imshow(Ref,[]);
    title("(c)")
    A = uint8(ImgSubstraction2(Img1,Img2,5,8,[512 512],@SSIAFA1));
    subplot(2,4,4);
    imshow(A,[]);
    title("(d)")
    A = uint8(ImgSubstraction2(Img1,Img2,5,8,[512 512],@SSIAFA2));
    subplot(2,4,5);
    imshow(A,[]);
    title("(e)")
    A = uint8(ImgSubstraction2(Img1,Img2,5,8,[512 512],@SSIAFA3));
    subplot(2,4,6);
    imshow(A,[]);
    title("(f)")
    A = uint8(ImgSubstraction2(Img1,Img2,5,8,[512 512],@SSIAFA4));
    subplot(2,4,7);
    imshow(A,[]);
    title("(g)")
end

function ShowOneAlgo(Img1,Img2,fun)
% Shows One Algorithms with Different Approximation Degrees
    Ref = imsubtract(Img1,Img2);
    subplot(2,4,1);
    imshow(Img1,[]);
    title("(a)")
    subplot(2,4,2);
    imshow(Img2,[]);
    title("(b)")
    subplot(2,4,3);
    imshow(Ref,[]);
    title("(c)")
    A = uint8(ImgSubstraction2(Img1,Img2,1,8,[512 512],fun));
    subplot(2,4,4);
    imshow(A,[]);
    title("(d)")
    A = uint8(ImgSubstraction2(Img1,Img2,2,8,[512 512],fun));
    subplot(2,4,5);
    imshow(A,[]);
    title("(e)")
    A = uint8(ImgSubstraction2(Img1,Img2,3,8,[512 512],fun));
    subplot(2,4,6);
    imshow(A,[]);
    title("(f)")
    A = uint8(ImgSubstraction2(Img1,Img2,4,8,[512 512],fun));
    subplot(2,4,7);
    imshow(A,[]);
    title("(g)")
    A = uint8(ImgSubstraction2(Img1,Img2,5,8,[512 512],fun));
    subplot(2,4,8);
    imshow(A,[]);
    title("(h)")
end

function [PSNR, MSSIM] = calculateMetrics(Img1,Img2,fun)
% Calculates PSNR and MSSIM for different Approximation Degrees for one
% Algorithm
    Ref = imsubtract(Img1,Img2);
    A = uint8(ImgSubstraction2(Img1,Img2,1,8,[512 512],fun));
    B = uint8(ImgSubstraction2(Img1,Img2,2,8,[512 512],fun));
    C = uint8(ImgSubstraction2(Img1,Img2,3,8,[512 512],fun));
    D = uint8(ImgSubstraction2(Img1,Img2,4,8,[512 512],fun));
    E = uint8(ImgSubstraction2(Img1,Img2,5,8,[512 512],fun));
    PSNR = [psnr(A,Ref);psnr(B,Ref);psnr(C,Ref);psnr(D,Ref);psnr(E,Ref)];
    MSSIM = [ssim(A,Ref);ssim(B,Ref);ssim(C,Ref);ssim(D,Ref);ssim(E,Ref)];
end

function Out = int2twoscomp(Value,n) 
% Returns Value in 2s Complement Binary, with n Bits in Vector Format
    BitString =  dec2bin(mod((Value),2^n),n);
    Out = uint8(flip(transpose(BitString-'0')));
end

function Output = ImgSubstraction2(Img1,Img2,k,n,size, fun)
% Uses ApprSubstraction for every Pixel in the two Grayscale Images
    Output = zeros(size(1),size(2));
    for i = (1:size(1))
        for j = (1:size(2))
            Output(i,j) = ApprSubstraction2(Img1(i,j),Img2(i,j),k,n, fun);
        end
    end
end

function Out = ApprSubstraction2(Int1,Int2,k,n, fun)
    % Calculates an n-Bit Substraction with k Approximated FA and (n-k) FA
    % k... Approximate Bits, n... Size
    Ain = int2twoscomp(Int1,n+1);
    Bin = int2twoscomp(-double(Int2),n+1);
    Cin = 0;
    Sum = zeros(n+1,1);
    if k>0
        for i = (1:k)
            [Cout, Sum(i)] =  fun(Ain(i),Bin(i),Cin);
            Cin = Cout;
        end
    end
    if k<8
        for j = (k+1:n+1)
            Sum(j) = xor(Ain(j),xor(Bin(j),Cin));
            Cin = (Ain(j) & Bin(j)) | (Cin & (Ain(j) | Bin(j)));
        end 
    end
    Out = bit2int(Sum,n+1,false);
    % Checks if Output would be negative and sets it to 0 if it is
    if (Out > Int1)
        Out = 0;
    end
end

function Output = ImgSubstraction(Img1,Img2,k,n,size, fun)
% Uses ApprSubstraction for every Pixel in the two Grayscale Images
    Output = zeros(size(1),size(2));
    for i = (1:size(1))
        for j = (1:size(2))
            Output(i,j) = ApprSubstraction(Img1(i,j),Img2(i,j),k,n, fun);
        end
    end
end

function Out = ApprSubstraction(Int1,Int2,k,n, fun)
    % Calculates an n-Bit Substraction with k Approximated FA and (n-k) FA
    % k... Approximate Bits, n... Size
    Ain = int2twoscomp(Int1,n);
    Bin = int2twoscomp(-double(Int2),n);
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
    Out = bit2int(Sum,n,false);
    % Checks if Output would be negative and sets it to 0 if it is
    if (Out > Int1)
        Out = 0;
    end
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