% s  = short (10:20)
% sm  = shortmedium (10:25)
% ss = very short (7.5:15)
% l  = long (10:40)
% ll = longdelayed (30:60)

index = 53
if index > numel(datasetSingleCells)
    datasetSingleCells{index}.CellID = '160715_AAAK';
    datasetSingleCells{index}.pos =  [908 340 86];
    datasetSingleCells{index}.posum =  [58 -183 78] ;
    datasetSingleCells{index}.VC70 = [1:3 11:13];
    datasetSingleCells{index}.VC70odor = {'sTrp' 'sArg' 'sAla' 'lTrp' 'lArg' 'lAla'};
    datasetSingleCells{index}.VC0 = [4:6 7:10] ;
    datasetSingleCells{index}.VC0odor = {'sTrp' 'sArg' 'sAla' 'lTrp' 'lArg' 'lAla' 'lAla'  };
    datasetSingleCells{index}.CCtest = [14 15];
    datasetSingleCells{index}.CCstim = [];
    datasetSingleCells{index}.CCodor = { };
    datasetSingleCells{index}.possiblyBadTrials = [2 9];
    datasetSingleCells{index}.comments = {'Strong and strange baseline fluctuation. Unclear whether biological or e.g. flow pulsing.'};
else
    disp('Really?')
end

[ numel(datasetSingleCells{index}.VC70), numel(datasetSingleCells{index}.VC70odor);
 numel(datasetSingleCells{index}.VC0), numel(datasetSingleCells{index}.VC0odor);
numel(datasetSingleCells{index}.CCstim), numel(datasetSingleCells{index}.CCodor)]

datasetSingleCells{index}.pos - datasetSingleCells{index-1}.pos

clear date
cd ..
save(strcat('dataset',date,'_',num2str(index),'.mat'),'datasetSingleCells')

