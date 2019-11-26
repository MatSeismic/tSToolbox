% save the results from linear inversion to txt file

load('LinearInv');

allX = xMat(:);
allY = yMat(:);
allT = MImodel(:);

[allLat, allLon] = minvtran(mstruct, allX, allY);
linearInv = [allLon allLat allT];

dlmwrite('WY.txt', linearInv);