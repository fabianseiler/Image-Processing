
%% FOM Calculation 
%Steps = [77,71,65,59];
Steps = [78, 73, 68, 63];

Energie_1 = 1.4509;
Energie_2 = 1.6694;
Energie_3 = 1.6678;
Energie_4 = 1.8697;
Energie_FA = 3.8435;
Energie_Extra = 0.805;

NMED1 = [0.00098,0.0044,0.0127,0.03];
NMED2 = [0.0039,0.0083,0.0155,0.0292];
NMED3 = [0.0039,0.0083,0.0155,0.0292];
NMED4 = [0.0039,0.0107,0.0237,0.0489];

for i = (1:4)
    Energie = Energie_4*i + Energie_Extra + Energie_FA*(8-i);
    FOM(i) = Steps(i)*Energie/(1-NMED4(i))
end

%% Quality Metrics for 8-Bit
for i = (1:4)
[MED(i) , NMED(i), MRED(i)] = MED_8Bit(5,8,@SSIAFA1);
end


%% Quality Metrics for 16/32-Bit
[MED, NMED, MRED] = MED_nBit(1,8,@SSIAFA1);


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

