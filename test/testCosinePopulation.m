% Tests for CosinePopulation

tauRef = .0025;
tauRC = .02;
n = 100;
intercepts = -1+2*rand(1,n);
maxRates = 100+200*rand(1,n);
V0 = 0;

sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
cp = CosinePopulation(1, sg, 'test');
do = addOrigin(cp, 'output', @(x) [x; x.^2]);

% addMerge(cp, [1 2]);
% assert(size(do.decoders, 2) == n+1)
% removeMerge(cp, 11);
% assert(size(do.decoders, 2) == n)

initClusters(cp)
initOriginErrors(cp);

setTolerance(do, [.1 .2].^2); 
tic, updateClusterIndices(cp); toc
cp.indices
error = getErrorEstimate(do, cp.indices)

% % this is not a fair comparison since error is wrt DIRECT ...
% f = zeros(1, length(cp.indices)); 
% for i = 1:length(f)
%     children = [];
%     if cp.indices(i) > n
%         children = DecodedOrigin.getChildren(cp.indices(i), cp.Z)
%     end
%     activeIndices = [cp.indices(i) setdiff(1:n, children)];
%     f(i) = plotDecodedOrigin(cp, 'output', [], 2, activeIndices); 
% end; 
% f


% indices = 2*n-1
% mse = [];
% for i = size(cp.Z, 1):-1:1
%     if mod(i, 5) == 0
%         mse = [mse plotDecodedOrigin(cp, 'output', [], 2, indices)];
%         pause
%         close(gcf)
%     end
%     indices = [setdiff(indices, n+i) cp.Z(i,:)]
% end
% figure, plot(mse)

% disp('Test OK')

% x = -1:.025:1;
% figure
% for i = n+1:2*n-1
%     i
%     hold off
%     merged = cp.Z(1:i-n,:);
%     merged = merged(:);
%     
%     plot(x, getRates(cp, x, 0, 0, setdiff(1:i, merged)), 'k');
%     hold on
% %     plot(x, getRates(cp, x, 0, 0, merged), 'k--');
%     plot(x, getRates(cp, x, 0, 0, cp.Z(i-n,:)), 'r');
%     plot(x, getRates(cp, x, 0, 0, i), 'g');
%     pause    
% end

% sg2 = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
% cp2 = CosinePopulation([1; 2], sg2, 'test');
% do2 = addOrigin(cp2, 'output', @(x) x(1,:) .* x(2,:));
% 
% initClusters(cp2)
% 
% 
% figure
% for i = n+1:2*n-1
%     i
%     hold off
%     merged = cp2.Z(1:i-n,:);
%     merged = merged(:);
%     cp2.encoders(:,merged)
%     
%     notmerged = setdiff(1:i, merged);
%     scatter(cp2.encoders(1,notmerged), cp2.encoders(2,notmerged), 'ko');
%     hold on
%     scatter(cp2.encoders(1,cp2.Z(i-n,:)), cp2.encoders(2,cp2.Z(i-n,:)), 'ro');
%     scatter(cp2.encoders(1,i), cp2.encoders(2,i), 'go');
%     pause    
% end
% 

