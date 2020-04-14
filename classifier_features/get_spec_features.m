
function f = get_spec_features(intCQT)

    vrs = (max(intCQT)./sum(intCQT));
    X = intCQT.^2;
    rms = sum(X);
    %FeatureSpectralCentroid
    vsc     = ([0:size(X,1)-1]*X)./(sum(X,1)+eps);
    
    %FeatureSpectralSpread    
    for (n = 1:size(X,2))
        vss(n)  = (([0:size(X,1)-1]-vsc(n)).^2*X(:,n))./sum(X(:,n));
    end
    vss     = sqrt(vss);

    X = intCQT;
    %SpectralSlope
    [vssl] = FeatureSpectralSlope (X, 0);
    %SpectralSkewness
    [vssk] = FeatureSpectralSkewness (X, 0);
    
       
    %SpectralRolloff
      % allocate memory
    vsr     = zeros(1,size(X,2));  
    %compute rolloff
    afSum   = sum(X,1);
    kappa   = 0.85;
    for (n = 1:length(vsr))
        vsr(n)  = find(cumsum(X(:,n)) >= kappa*afSum(n), 1); 
    end
    
    [vsf] = FeatureSpectralFlux (X, 0);
    [vtf] = FeatureSpectralFlatness (X, 0);
    [vsd] = FeatureSpectralDecrease (X, 0);
    [vtsc] = FeatureSpectralCrestFactor (X, 0);
    
    f = [vrs' rms' vsc' vss' vssl' vssk' vsr' vsf' vtf' vsd' vtsc']; 
end
% vrs = ((max(intCQT)./sum(intCQT))>.1);
% vrs = ((max(intCQT)./sum(intCQT))>.02);
% sumx = diag(vrs);
