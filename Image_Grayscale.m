%% Grayscale Conversion
rgb_img = uint8(imread("peppers.png"));


%ShowAllAlgos(rgb_img,3)
%ShowOneAlgo(rgb_img,@SSIAFA4)
%[PSNR, MSSIM] = calculateMetrics(rgb_img,@SSIAFA3)

%% functions
function [PSNR, MSSIM] = calculateMetrics(rgbImg,fun)
% Calculates PSNR and MSSIM for different Approximation Degrees for one
% Algorithm
    Ref = uint8(GrayscaleConversion(rgbImg,0,8,fun));
    A = uint8(GrayscaleConversion(rgbImg,1,8,fun));
    B = uint8(GrayscaleConversion(rgbImg,2,8,fun));
    C = uint8(GrayscaleConversion(rgbImg,3,8,fun));
    D = uint8(GrayscaleConversion(rgbImg,4,8,fun));
    E = uint8(GrayscaleConversion(rgbImg,5,8,fun));
    PSNR = [psnr(A,Ref);psnr(B,Ref);psnr(C,Ref);psnr(D,Ref);psnr(E,Ref)];
    MSSIM = [ssim(A,Ref);ssim(B,Ref);ssim(C,Ref);ssim(D,Ref);ssim(E,Ref)];
end

function ShowOneAlgo(rgbImg,fun)
% Shows One Algorithms with Different Approximation Degrees
    subplot(3,3,1)
    imshow(rgbImg,[])
    title("(a)")
    Ref = GrayscaleConversion(rgbImg,0,8,fun);
    subplot(3,3,2)
    imshow(Ref,[])
    title("(b)")
    A = GrayscaleConversion(rgbImg,1,8,fun);
    subplot(3,3,3)
    imshow(A,[])
    title("(c)")
    A = GrayscaleConversion(rgbImg,2,8,fun);
    subplot(3,3,4)
    imshow(A,[])
    title("(d)")
    A = GrayscaleConversion(rgbImg,3,8,fun);
    subplot(3,3,5)
    imshow(A,[])
    title("(e)")
    A = GrayscaleConversion(rgbImg,4,8,fun);
    subplot(3,3,6)
    imshow(A,[])
    title("(f)")
    A = GrayscaleConversion(rgbImg,5,8,fun);
    subplot(3,3,7)
    imshow(A,[])
    title("(g)")
end

function ShowAllAlgos(rgbImg,k)
% Shows all the Different Image Addition Results with all Algorithms with
% an Approximation Degree of k/8 
    subplot(2,3,1)
    imshow(rgbImg,[])
    title("(a)")
    Ref = GrayscaleConversion(rgbImg,0,8,@SSIAFA1);
    subplot(2,3,2)
    imshow(Ref,[])
    title("(b)")
    A = GrayscaleConversion(rgbImg,k,8,@SSIAFA1);
    subplot(2,3,3)
    imshow(A,[])
    title("(c)")
    A = GrayscaleConversion(rgbImg,k,8,@SSIAFA2);
    subplot(2,3,4)
    imshow(A,[])
    title("(d)")
    A = GrayscaleConversion(rgbImg,k,8,@SSIAFA3);
    subplot(2,3,5)
    imshow(A,[])
    title("(e)")
    A = GrayscaleConversion(rgbImg,k,8,@SSIAFA4);
    subplot(2,3,6)
    imshow(A,[])
    title("(f)")
end

function Output = GrayscaleConversion(rgbImg,k,n,fun)
% Converts an RGB Image to a Grayscale Image with an RCA that has k
% Approximated FA and (n-k) Exakt FA
    RedValues = uint8(rgbImg(:,:,1))/3;
    GreenValues = uint8(rgbImg(:,:,2))/3;
    BlueValues = uint8(rgbImg(:,:,3))/3;
    ImgSize = size(rgbImg);
    RGValues = ImgAddition(RedValues,GreenValues,k,n,ImgSize,fun);
    Output = ImgAddition(RGValues,BlueValues,k,n,ImgSize,fun);
end

function Output = ImgAddition(Img1,Img2,k,n,size,fun)
% Uses ApprAddition for every Pixel for two different Image Matricies
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
    Out = bit2int(Sum,n,false);
    % Checks the Edge Case where the Rounded Values Reach over 2^n-1 and
    % sets the value back to 2^n-1
    if Out > 2^n-1
        Out = 2^n-1;
    end
    Out = uint8(Out);
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