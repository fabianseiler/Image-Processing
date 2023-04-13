%% Description

% Old file that isnt used anymore
% Includes the first trys of Image Processing


%% Image Addition
Img1 = imread('rice.png');
Img2 = imread('cameraman.tif');
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


%title(['Appr/Exact: ', num2str(i),'/',num2str((8-i)),]);
%[MSSIM, ~] = ssim(A,Ref);
%xlabel(['PSNR: ',num2str(psnr(A,Ref)), ', SSIM: ',num2str(2), ', MSSIM: ', num2str(MSSIM)])


%% Image Addition mit Hardcoded RCA
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
[MSSIM, SSIMmap] = ssim(A,Ref);
xlabel(['PSNR: ', num2str(psnr(A,Ref)), ', MSSIM: ' , num2str(MSSIM)])
figure
imshow(Ref,[])



%% Image Substraction
Img1 = imread('motion01.512.tiff');
subplot(2,2,1);
imshow(Img1,[])
Img2 = imread('motion05.512.tiff');
subplot(2,2,2);
imshow(Img2,[])
Ref = imsubtract(uint8(Img2)/2,uint8(Img1)/2);
subplot(2,2,3);
imshow(Ref,[])
A = uint8(ImgSubstraction(uint8(Img1)/2,uint8(Img2)/2,0,8,[512 512],@SSIAFA1));
subplot(2,2,4);
imshow(A,[])

%% Grayscale Conversion
% Genaue Berechnung anschauen
PSNR = zeros(6,1);
SSIM = zeros(6,1);
rgb_img = imread('peppers.png');
size(rgb_img)
%rgb_img = imresize(rgb_img,[256 256]);
Ref = uint8(ImgAddition(uint8(rgb_img(:,:,1))/3,ImgAddition(uint8(rgb_img(:,:,2))/3,uint8(rgb_img(:,:,3))/3,0,8,[384 512],@SSIAFA3),0,8,[384 512],@SSIAFA3));
subplot(2,4,1);
imshow(rgb_img,[])
subplot(2,4,2);
imshow(Ref,[])
for i = (0:5)
    A = uint8(ImgAddition(uint8(rgb_img(:,:,1))/3,ImgAddition(uint8(rgb_img(:,:,2))/3,uint8(rgb_img(:,:,3))/3,i,8,[384 512],@SSIAFA3),i,8,[384 512],@SSIAFA3));
    subplot(2,4,i+3);
    imshow(A,[]);
    title(['Appr/Exact: ', num2str(i),'/',num2str((8-i)),]);
    [MSSIM, SSIMmap] = ssim(A,Ref);
    xlabel(['PSNR: ',num2str(psnr(A,Ref)), ', SSIM: ',num2str(2), ', MSSIM: ', num2str(MSSIM)]);
end

%% Image Multiplication
Img1 = imread('cameraman.tif');
Img2 = imread('pout.tif');
Img2 = imresize(Img2, [256 256]);
A = immultiply(Img1/8,Img2/8);
imshow(A,[])

%% Quality Metrics for 8-Bit
for i = (1:1)
[MED(i) , NMED(i), MRED(i)] = MED_8Bit(5,8,@SSIAFA1);
end



%% Quality Metrics for 16/32-Bit
[MED, NMED, MRED] = MED_nBit(1,8,@SSIAFA1);


%% FOM Calculation 
%Steps = [77,71,65,59];
Steps = [78, 73, 68, 63];

Energie_1 = 1.4509;
Energie_2 = 1.6694;
Energie_3 = 1.6678;
Energie_4 = 1.8697;
Energie_FA = 9.87;
Energie_Extra = 1.33;

NMED1 = [0.00098,0.0044,0.0127,0.0313];
NMED2 = [0.0039,0.0083,0.0155,0.0292];
NMED3 = [0.0039,0.0083,0.0155,0.0292];
NMED4 = [0.0039,0.0107,0.0237,0.0489];

for i = (1:4)
    Energie = Energie_4*i + 1.33 + 9.87*(8-i); 
    FOM(i) = Steps(i)*Energie/(1-NMED4(i));
end


%% Testing Area
flip(int2twoscomp(55,8))


%% functions
function [MED, NMED, MRED] = MED_8Bit(k,n,fun)
    ED_Sum = 0;
    for a = (0:1:2^n-1)
        for b = (0:1:2^n-1)
            Sum_exact = ApprAddition(a,b,0,n, fun);
            Sum_appr = ApprAddition(a,b,k,n, fun);
            ED_Sum = ED_Sum + abs(double(Sum_exact)-double(Sum_appr));
            if(Sum_exact == 0)
                RED(a+1,b+1) = 0;
            else
                RED(a+1,b+1) = double(abs(double(Sum_exact)-double(Sum_appr))/double(Sum_exact)); 
            end
        end
    end
    MED = ED_Sum/(2^(2*n));
    NMED = MED/(2^n-1);
    MRED = sum(RED,"all")/(2^(2*n));
end

function [MED, NMED, MRED] = MED_nBit(k,n,fun)
% Hat einen Error bei der Berechnung
    InputValues1 = randi(2^n-1,[1000,1]);
    InputValues2 = randi(2^n-1,[1000,1]);
    ED_Sum = 0;
    RED = 0;
    for a = InputValues1(:)'
        for b = InputValues2(:)'
            Sum_exact = ApprAddition(a,b,k,n, fun);
            Sum_appr = ApprAddition(a,b,k,n, fun);
            ED_Sum = ED_Sum + abs(double(Sum_exact)-double(Sum_appr));
            if(Sum_exact ~= 0)
               RED = RED + double(abs(double(Sum_exact)-double(Sum_appr))/double(Sum_exact));
            end   
        end
    end
    MED = ED_Sum/(1000^2-1);
    NMED = MED/(2^n-1);
    MRED = RED/(1000^2-1);
end

function Out = int2twoscomp(Value,n) 
% Returns -Value in 2s Complement Binary, with n+1 Bits

%MSB shiften und auf 1 setzen
     BitString =  dec2bin(mod((-Value),2^n),n);
     Out = transpose(BitString-'0');
end

function Output = ImgSubstraction(Img1,Img2,k,n,size, fun)
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
    Ain = int2bit(Int1,n,false);
    Bin = uint8(flip(int2twoscomp(Int2,n)));
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

function Out= RCA(Int1, Int2, fun)
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


