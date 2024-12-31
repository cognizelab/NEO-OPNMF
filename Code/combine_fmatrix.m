function [Model,EVA,Log,Data] = combine_fmatrix(result1,result2,output)
%% Adding more CVs to enhance the stability.
load(result2);
EVA1 = EVA; Log1 = Log;
load(result1);
%%
EVA.increased_SED = [EVA.increased_SED;EVA1.increased_SED];
EVA.increased_RMSE = [EVA.increased_RMSE;EVA1.increased_RMSE];
EVA.VI = [EVA.VI;EVA1.VI];
EVA.CI = [EVA.CI;EVA1.CI];
EVA.MV = cat(3,EVA.MV,EVA1.MV);
EVA.Jacard = [EVA.Jacard;EVA1.Jacard];
EVA.Rand = [EVA.CI;EVA1.Rand];
EVA.zRand = [EVA.zRand;EVA1.zRand];
EVA.aRand = [EVA.aRand;EVA1.aRand];

Log.CV.n = Log.CV.n + Log1.CV.n;

if nargin < 3 || isempty(output)
    Log.OS.filename = 'FMATRIX';
    save FMATRIX Model EVA Log Data
else
    Log.OS.filename = output;
    eval(['save ',output,' Model EVA Log Data']);
end

end