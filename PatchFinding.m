function Ans = PatchFinding(Images, indices, Fsel,patchSize,refInt, searchArea, neighbourhood, type)
    [H,W,C,K] = size(Images);
    Ans= zeros(H,W,C,K);
    Count= zeros(H,W,C,K);
    N= patchSize;                %size of ref patch
    G= refInt;                %samlping interval
    for i=1:K
        for ypos=1:G:H-N+1
            for xpos=1:G:W-N+1
                i
                ypos
                xpos
                RefP= Images(ypos: ypos+N-1, xpos: xpos+N-1, :, i);
                RefPos= [i, ypos, xpos];
                [Pjk, Pos]= Match(Images, RefP , Fsel, RefPos, searchArea, neighbourhood, type);
                
                Omega= selectReliable(Pjk, indices, Pos, patchSize);
                Q= Denoise(Pjk, Omega);
                
                for k=1: size(Q, 2)
                    patch= reshape(Q(:,k), N,N, C);
                    Ans(Pos(k,2):Pos(k,2)+N-1, Pos(k,3): Pos(k,3)+N-1, :, Pos(k,1))= Ans(Pos(k,2):Pos(k,2)+N-1, Pos(k,3): Pos(k,3)+N-1, :, Pos(k,1))+ patch;
                    Count(Pos(k,2):Pos(k,2)+N-1, Pos(k,3): Pos(k,3)+N-1, :, Pos(k,1))= Count(Pos(k,2):Pos(k,2)+N-1, Pos(k,3): Pos(k,3)+N-1, :, Pos(k,1))+1;
                end                
            end
        end
    end
    Count= max(Count,1);
    Ans= Ans./Count;
end

function [Pjk, Pos] = Match(Images, RefP, Fsel,RefPos, searchArea, neighbourhood, type)
    [~,~,~,K] = size(Images);
    %Pos are the positions (Frame no, x, y) of those matched patches
    Ind=1;
    for i= max([1,RefPos(1)-neighbourhood]):min([K,RefPos(1)+neighbourhood])
        [P, Posf]= Select(Images(:,:, :, i), RefP, Fsel,RefPos, searchArea, type);
        Pjk(:, Ind: Ind+Fsel-1)= P;
        Pos(Ind: Ind+Fsel-1, 2:3)= Posf;
        Pos(Ind: Ind+Fsel-1, 1)= i;
        Ind= Ind+ Fsel;
    end
end

function [P, Posf] = Select(Image, RefP, Fsel,RefPos, searchArea, type)
    [H,W,C]=size(Image);
    N= size(RefP, 1);
    P= zeros(N*N*C, Fsel);
    Allmads= zeros(0);
    Allpos= zeros(0);
    if strcmp(type,'exhaustive')
        for i=max([1,RefPos(2)-searchArea]): min([H-N+1,RefPos(2)+searchArea])
            for j=max([1,RefPos(3)-searchArea]): min([W-N+1,RefPos(3)+searchArea])
                currPatch= Image(i: i+N-1, j: j+N-1, :);
                MAD= abs(RefP-currPatch);
                Allmads(end+1)= sum(MAD,'all');
                Allpos(end+1, :)= [i j];
            end
        end
    elseif strcmp(type,'fast')
        assert(mod(N,2) == 0,'Ref Patch is odd');
        assert(mod(H,2) == 0 && mod(W,2) == 0,'Frame size is odd');
        RefP_labels = repmat(cast([1 2; 3 4],'uint8'),N/2,N/2,C);
        Image_labels = repmat(cast([1 4; 2 3],'uint8'),H/2,W/2,C);
        Allmads_temp = [];
        Allpos_temp = [];
        for i=max([1,RefPos(2)-searchArea]): min([H-N+1,RefPos(2)+searchArea])
            for j=max([1,RefPos(3)-searchArea]): min([W-N+1,RefPos(3)+searchArea])
                currPatch= Image(i: i+N-1, j: j+N-1, :);
                Patch_label = Image_labels(i,j,:);
                Pattern = (RefP_labels == Patch_label);
                MAD = abs(RefP(Pattern)-currPatch(Pattern));
                X = inf(1,4);
                X(Patch_label) = sum(MAD,'all');
                Allmads_temp = [Allmads_temp; X];
                Allpos_temp = [Allpos_temp; [i j]];
            end
        end
        [~,Inds_temp] = mink(Allmads_temp,Fsel);
        for i = 1:4
            Pos_t = Allpos_temp(Inds_temp(:,i),:);
            for j = 1:Fsel
                currPatch = Image(Pos_t(j,1):Pos_t(j,1)+N-1,Pos_t(j,2):Pos_t(j,2)+N-1);
                MAD = sum(abs(currPatch-RefP),'all');
                Allmads(end+1) = MAD;
                Allpos(end+1,:) = Pos_t(j,:);
            end
        end
    else
        throw('Patch Finding type not found')
    end
    [~, Inds]= mink(Allmads, Fsel);
    Posf= Allpos(Inds, :);
    for i=1: Fsel
        currPatch= Image(Posf(i,1): Posf(i,1)+N-1, Posf(i,2): Posf(i,2)+N-1, :);
        P(:,i) = reshape(currPatch, N*N*C, 1);
    end
    assert(~any(isnan(P),'all'),"Output array contains inf values, please use bigger frames");
end

function Omega= selectReliable(Pjk, indices, Pos, patchSize)
    [N,M] = size(Pjk);
    S= mean(Pjk, 2);
    temp = std(Pjk,0,2);
    temp= temp .^ 2;
    sigma = sum(temp, 'all')/N;
    sigma = 2* sqrt(sigma);
    Omega= abs(Pjk-S)<= sigma;
    Omega2= zeros(N,M);
    for k=1: size(Pjk, 2)
        patch=indices(Pos(k,2):Pos(k,2)+patchSize-1, Pos(k,3): Pos(k,3)+patchSize-1, :, Pos(k,1));
        Omega2(:,k)= patch(:);
    end
    Omega= Omega2 & Omega;    
end
