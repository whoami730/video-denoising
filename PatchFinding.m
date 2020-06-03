function Ans = PatchFinding(Images,Fsel,patchSize,refInt,type)
    [H,W,C,K] = size(Images);
    Ans= zeros(H,W,C,K);
    Count= zeros(H,W,C,K,'uint8');
    N= patchSize;                %size of ref patch
    G= refInt;                %samlping interval
    for i=1:K
        for ypos=1:G:H-N+1
            for xpos=1:G:W-N+1
                for ch=1:C
                    RefP= Images(ypos: ypos+N-1, xpos: xpos+N-1, ch, i);
                    [Pjk, Pos]= Match(Images, RefP, ch, Fsel,type);
                    Omega= selectReliable(Pjk);
                end
            end
        end
    end
end

function [Pjk, Pos] = Match(Images, RefP, ch, Fsel,type)
    [~,~,~,K] = size(Images);
    N= size(RefP, 1);
    Pjk= zeros(N*N, Fsel*K);    %Pjk is the matrix of matched patches in the chanel ch
    Pos= zeros(Fsel*K, 3,'uint8');      %Pos are the positions (Frame no, x, y) of those matched patches
    for i= 1:K
        [P, Posf]= Select(Images(:,:, ch, i), RefP, Fsel,type);
        Pjk(:, (i-1)*Fsel+1: i*Fsel)= P;
        Pos((i-1)*Fsel+1: i*Fsel, 2:3)= Posf;
        Pos((i-1)*Fsel+1: i*Fsel, 1)= i;
    end
end

function [P, Posf] = Select(Image, RefP, Fsel,type)
    [H,W]=size(Image);
    N= size(RefP, 1);
    P= zeros(N*N, Fsel);
    Ar= inf(Fsel, 3);
    Posf= zeros(Fsel, 2,'uint8');
    if strcmp(type,'exhaustive')
        for i=1: H-N+1
            for j=1:W-N+1
                currPatch= Image(i: i+N-1, j: j+N-1);
                MAD= abs(RefP-currPatch);
                currVal= [sum(MAD, 'all') i j];
                for k=1:Fsel
                    if Ar(k,1) > currVal(1)
                        temp= Ar(k,:);
                        Ar(k,:)= currVal;
                        currVal= temp;
                    end
                end
            end
        end
    elseif strcmp(type,'fast')
        assert(mod(N,2) == 0,'Ref Patch is odd');
        assert(mod(H,2) == 0 && mod(W,2) == 0,'Frame size is odd');
        RefP_labels = repmat(cast([1 2; 3 4],'uint8'),N/2);
        Image_labels = repmat(cast([1 4; 2 3],'uint8'),H/2,W/2);
        Ar_temp = inf(Fsel*4,3);
        for i = 1:H-N+1
            for j = 1:W-N+1
                currPatch= Image(i: i+N-1, j: j+N-1);
                Patch_label = Image_labels(i,j);
                Pattern = (RefP_labels == Patch_label);
                MAD = abs(RefP(Pattern)-currPatch(Pattern));
                currVal= [sum(MAD, 'all') i j];
                for k=1:Fsel*4
                    if Ar_temp(k,1) > currVal(1)
                        temp= Ar_temp(k,:);
                        Ar_temp(k,:)= currVal;
                        currVal= temp;
                    end
                end
            end
        end
        for i = 1:Fsel*4
            currVal = Ar_temp(i,:);
            a = currVal(2); b = currVal(3);
            currVal(1) = sum(abs(Image(a:a+N-1,b:b+N-1) - RefP),'all');
            for j = 1:Fsel
                if Ar(j,1) > currVal(1)
                    temp = Ar(j,:);
                    Ar(j,:) = currVal;
                    currVal = temp;
                end
            end
        end
    else
        throw('Patch Finding type not found')
    end
    for i=1: Fsel
        Posf(i,:) = Ar(i,2:3);
        currPatch= Image(Ar(i,2): Ar(i,2)+N-1, Ar(i,3): Ar(i,3)+N-1);
        P(:,i) = reshape(currPatch, N*N, 1);
    end
    assert(~any(isnan(P),'all'),"Output array contains inf values, please use bigger frames");
end

function Omega= selectReliable(Pjk)
    [N,M] = size(Pjk);
    Pjk= double(Pjk);
    S= sum(Pjk, 2)/M;
    D= Pjk-S;
    D= D.*D;
    Ds=sum(D,2)/(M-1);
    sigma= 2*sqrt(sum(Ds, 'all'));
    Omega= zeros(N,M);
    for i=1: N
        for j=1: M
            if (abs(Pjk(i,j)-S(i))<= sigma)
                Omega(i,j)=1;
            end
        end
    end
end
