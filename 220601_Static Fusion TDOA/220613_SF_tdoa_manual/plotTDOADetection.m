function plotTDOADetection(tdoaDetections)
    x = zeros(0,1);
    y = zeros(0,1);
    for i = 1:numel(tdoaDetections)
        [thisX, thisY] = tdoaHyperbola(tdoaDetections{i},'2D');
        x = [x;nan;thisX(:)]; %#ok<AGROW>
        y = [y;nan;thisY(:)]; %#ok<AGROW>
    end

    linkdata
    plot(x,y)
    xlim([-10e3 10e3]); ylim([-10e3 10e3]);
    drawnow
    
end