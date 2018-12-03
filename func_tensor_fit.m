function [] = func_tensor_fit(inputArg1)
% Create functional correlational tensors using fMRI data. Takes timeseries data extracted using 3dNetCorr; alternatively, any form of
% [X1,Y1,Z1, t1,t2,t3,t4...] should work.

%   Create Functional Correlational Tensors using fMRI data:

%load the data
ds = datastore(inputArg1); %name of file goes here
tt = tall(ds);
raw = table2array(gather(tt));
voxels = raw(:,1:3);
data = raw(:,4:153);

%create the design matrix, Mi
M = nan(26,3);
Mi = nan(26,6);
Morig = nan(26,3);
i = 1;
unit = [1 0 -1];
for c = unit(1:3)
    for b = unit(1:3)
        for a = unit(1:3)
            M(i,:) = [a b c];
            Morig(i,:) = [a b c];
            if sum(abs(M(i,:)))' == 2
                M(i,:) = M(i,:).*[(1/sqrt(2)) (1/sqrt(2)) (1/sqrt(2))];
            elseif sum(abs(M(i,:)))' == 3
                M(i,:) = M(i,:).*[(1/sqrt(3)) (1/sqrt(3)) (1/sqrt(3))];
            end
            Mi(i,:) = [M(i,1)^2 2*M(i,1)*M(i,2) 2*M(i,1)*M(i,3) M(i,2)^2 2*M(i,2)*M(i,3) M(i,3)^2]; %normalize
            i = i+1;
            [a;b;c];
        end
    end
end

Mi(14,:) = []
Morig(14,:) = []

%beefjj = 
results = zeros(length(voxels),6);
for jj = 1:length(voxels)
    %if voxels(jj,1) > min(voxels(:,1)) & voxels(jj,2) > min(voxels(:,2)) & voxels(jj,3) > min(voxels(:,3)) &...
    %  voxels(jj,1) < max(voxels(:,1)) & voxels(jj,2) < max(voxels(:,2)) & voxels(jj,3) < max(voxels(:,3))
    C = nan(26,1);
    for ii = 1:26
        try
            add1 = Morig(ii,1);
            add2 = Morig(ii,2);
            add3 = Morig(ii,3);
            R = corrcoef(data(jj, :),...
                data(voxels(:,1) == (voxels(jj,1)+add1) & ...
                voxels(:,2) == (voxels(jj,2)+add2) & ...
                voxels(:,3) == (voxels(jj,3)+add3),:));
            C(ii) = R(2)^2;
        catch
            C(ii) = 0;
        end
    end
    
    T = inv(Mi' * Mi) * Mi' * C;
    results(jj,:) = T';
    
    
end
results_1 = [voxels results(:,1)];
results_2 = [voxels results(:,2)];
results_3 = [voxels results(:,3)];
results_4 = [voxels results(:,4)];
results_5 = [voxels results(:,5)];
results_6 = [voxels results(:,6)];
tensorname = '_tensor'
save(strcat(inputArg1,tensorname,'_1'),'results_1','-ascii','-double')
save(strcat(inputArg1,tensorname,'_2'),'results_2','-ascii','-double')
save(strcat(inputArg1,tensorname,'_3'),'results_3','-ascii','-double')
save(strcat(inputArg1,tensorname,'_4'),'results_4','-ascii','-double')
save(strcat(inputArg1,tensorname,'_5'),'results_5','-ascii','-double')
save(strcat(inputArg1,tensorname,'_6'),'results_6','-ascii','-double')



end
