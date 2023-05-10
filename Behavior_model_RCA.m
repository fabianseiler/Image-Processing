function [Cout, Sum_n] = ApprAddition(Int1, Int2, Cin, k, n, fun)
    % Calculates an n-Bit Addition with k Approximated FA and (n-k) exact FA
    % Int1, Int2 are the two input numbers in range (0:2^n - 1) 
    % Cin ... Bit Value of Carry In in range (0:1)
    % k... Approximate Bits, n... Size
    % fun ... Input function for approximated algorithm
    % Cout ... Carry Out Bit
    % Sum_n ... Sum of the Calculation with n bits
    Ain = int2bit(Int1,n,false);
    Bin = int2bit(Int2,n,false);
    Sum = zeros(n,1);
    if k>0
        for i = (1:k)
            % Approximated Algorithm
            [Cout, Sum(i)] =  fun(Ain(i),Bin(i),Cin);
            Cin = Cout;
        end
    end
    if k<n
        for j = (k+1:n)
            % Exact Full Adder logic functions
            Sum(j) = xor(Ain(j),xor(Bin(j),Cin));
            Cin = (Ain(j) & Bin(j)) | (Cin & (Ain(j) | Bin(j)));
        end 
    end
    Sum_n = bit2int(Sum,n,false);
    Cout = Cin;
end

function [Cout, Sum] = Algorithm1(Ain,Bin,Cin)
    Cout= (Ain&Bin) | Cin;
    Sum = ~Cout;
end

function [Cout, Sum] = Algorithm2(Ain,Bin,Cin)
    Cout= (Ain&Cin) | Bin;
    Sum = ~Cout;
end

function [Cout, Sum] = Algorithm3(Ain,Bin,Cin)
    Cout= (Bin&Cin) | Ain;
    Sum = ~Cout;
end

function [Cout, Sum] = Algorithm4(Ain,Bin,Cin)
    Cout= (Ain|Bin)&Cin;
    Sum = ~Cout;
end
