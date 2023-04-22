%% FOM Calculation 
%Steps = [77, 71, 65, 59, 53];         %Steps for Algo1
Steps = [78, 73, 68, 63, 58];       %Steps for Algo2,3,4

%Calculated Energies
Energie_1 = 1.4509;
Energie_2 = 1.6694;
Energie_3 = 1.6678;
Energie_4 = 1.8697;
Energie_FA = 3.8435;
Energie_Extra = 0.805;

%Calculated NMEDs
NMED1 = [0.00098,0.0044,0.0127,0.03,0.0702];
NMED2 = [0.0039,0.0083,0.0155,0.0292,0.0556];
NMED3 = [0.0039,0.0083,0.0155,0.0292,0.0556];
NMED4 = [0.0039,0.0107,0.0237,0.0489,0.0973];

FOM = zeros(5,1);
for i = (1:5)
    Energie = Energie_4*i + Energie_Extra + Energie_FA*(8-i);
    FOM(i) = Steps(i)*Energie/(1-NMED4(i));
end


%% Average Quality Metrics
PSNR_Add = [[51.1243,51.1243,51.1243,51.1243];
            [48.4595,47.1640,47.1695,48.1774];
            [43.7922,43.3313,43.3060,43.4088];
            [38.0888,38.2048,38.3122,37.6812];
            [32.0605,32.9769,32.9275,32.0572]];
MSSIM_Add = [[0.9976,0.9976,0.9976,0.9976];
             [0.9957,0.9940,0.9941,0.9952];
             [0.9884,0.9858,0.9860,0.9865];
             [0.9619,0.9583,0.9603,0.9576];
             [0.8966,0.8901,0.8934,0.8920]
            ];
PSNR_Sub = [[54.2460,52.2109,52.2117,53.5533];
            [45.8122,45.7904,45.8030,45.7522];
            [40.4358,40.6885,40.5615,40.8920];
            [35.3346,35.6264,35.7504,36.0118];
            [30.5750,30.7005,30.9533,31.3007]
           ];
MSSIM_Sub = [[0.9875,0.9755,0.9757,0.9879];
             [0.9157,0.9251,0.9250,0.9350];
             [0.8157,0.8198,0.8181,0.8762];
             [0.7124,0.6977,0.6996,0.8097];
             [0.6346,0.6043,0.6026,0.7849]
            ];
PSNR_Gray = [[57.1414,53.0294,52.7967,54.1623];
             [51.3954,49.5837,49.9157,48.9920];
             [44.9026,45.1002,45.6004,42.9854];
             [37.6214,40.0990,39.9520,36.5704];
             [29.7620,34.0016,33.3269,30.4971]
            ];
MSSIM_Gray = [[0.9989,0.9977,0.9977,0.9977];
              [0.9961,0.9957,0.9959,0.9943];
              [0.9841,0.9871,0.9884,0.9818];
              [0.9394,0.9609,0.9638,0.9472];
              [0.8168,0.8899,0.8986,0.8627]
             ];

PSNR_AVG = (PSNR_Add + PSNR_Sub + PSNR_Gray)/3;
MSSIM_AVG = (MSSIM_Add + MSSIM_Sub + MSSIM_Gray)/3;

%% Quality Metrics for 8-Bit
for i = (1:4)
[MED(i) , NMED(i), MRED(i)] = MED_8Bit(5,8,@SSIAFA1);
end


%% Quality Metrics for 16/32-Bit
[MED, NMED, MRED] = MED_nBit(1,8,@SSIAFA1);


%% functions
function [MED, NMED, MRED] = MED_8Bit(k,n,fun)
% Calculates Error Metrics for an 8Bit RCA with k ApprFa and (n-k) Exakt FA
% with the function fun for the Logic Outputs of the Approximated FA
% Error Metrics are Calculated as in the Literature
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
% Creates an random 1000x1000 Input Array and Calulates the Error Metrics
% for n Bits. 
% Uses the Logic Function of ApprFA as fun
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
    % Reference ApprFA Algorithms

    
    %-------SIAFA1------------
    %Sum = ~Bin | (~Ain & ~Cin);
    %Cout = ~Sum;

    %-------SIAFA4------------
    Sum = ~Cin | (~Ain & ~Bin);
    Cout = ~Sum; 
    %%
end

