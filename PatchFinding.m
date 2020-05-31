function Ans = PatchFinding(Images)
K= size(Images, 4);  %number of frames 
H= size(Images, 1);  %Height
W= size(Images, 2);  %Width
C= size(Images, 3);  %Channels
Fsel=5;              %number of matching patches selected per frame for ref patch
Ans= double(zeros(H,W,C,K));
Count= zeros(H,W,C,K);
N= 8;                %size of ref patch
G= 4;                %samlping interval
for i=1:K
    for ypos=1:4:H-G+1
        for xpos=1:4:W-G+1
            for ch=1:C
                RefP= Images(ypos: ypos+N-1, xpos: xpos+N-1, ch, i);
                [Pjk, Pos]= Match(Images, RefP, ch, Fsel);
            end
            
        end
    end
end
end

function [Pjk, Pos] = Match(Images, RefP, ch, Fsel)
K= size(Images, 4);  %number of frames 
H= size(Images, 1);  %Height
W= size(Images, 2);  %Width
C= size(Images, 3);  %Channels
N= size(RefP, 1);
Pjk= zeros(N*N, Fsel*K);    %Pjk is the matrix of matched patches in the chanel ch
Pos= zeros(Fsel*K, 3);      %Pos are the positions (Frame no, x, y) of those matched patches
for i= 1:K
    [P, Posf]= Select(Images(:,:, ch, i), Refp, Fsel);
    Pjk(:, i*(Fsel-1)+1: i*Fsel)= P;
    Pos(i*(Fsel-1)+1: i*Fsel, 2:3)= Posf;
    Pos(i*(Fsel-1)+1: i*Fsel, 1)= K;
end
end

function [P, Posf] = Select(Image, RefP, Fsel)
H=size(Image, 1);
W=size(Image, 2);
N= size(RefP, 1);
P= zeros(N*N, Fsel);
Ar= N*N*1000*ones(Fsel, 3);
Posf= zeros(Fsel, 3);
for i=1: H-N+1
    for j=1:W-N+1
        currPatch= Image(i: i+N-1, j: j+N-1);
        MAD= RefP-currPatch;
        MAD= abs(MAD);
        currVal= zeros(1,3);
        currVal(1)= sum(MAD, 'all');
        currVal(2)= i;
        currVal(3)= j;
        for k=1:Fsel
            if Ar(k,1) > currVal(1)
                temp= Ar(k);
                Ar(k)= currVal;
                currVal= temp;
                
            end
        end
    end
end
for i=1: Fsel
    Posf(i,1)= Ar(i,2);
    Posf(i,2)= Ar(i,3);
    currPatch= Image(Ar(i,2): Ar(i,2)+N-1, Ar(i,3): Ar(i,3)+N-1);
    curr= reshape(currPatch, N*N, 1);
    P(i)= curr;
end
end


