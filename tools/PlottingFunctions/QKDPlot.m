classdef QKDPlot
    % QKDPlot Helper class that defines methods to assist with plotting.
    %
    % If the selected parameter to plot is not given by scalar real values,
    % the plotting function instead will replace it with the index of the
    % data point. If multiple results are being plotted at once, then each
    % set of results must have the same number of data points.
    % Each axis can use 1 of 3 styles:
    % * linear: The axis is a good old fashion linear axis.
    % * log: The axis is log scaled.
    % * dB: The data is converted to dB using -10*log10(x). Note this only
    %   works when the data is scalar real values. (non positive values
    %   will be ignored for plotting).
    properties (Constant)
        keyRateTag = "_keyRate_"; % key value to indicate that you want to plot the key rate and not a parameter.
    end

    methods(Static)
        function simple1DPlot(qkdInput,results,options)
            % simple1DPlot Auto plot a single scan paramter vs key rate.
            %
            % Input:
            % * qkdInput: QKDSolverInput used to calculate the results. It
            %   must have exactly 1 scan parameter.
            % * results: The corresponding results 1D array for the
            %   qkdInput.
            % * xScaleStyle ("linear"): Name-value argument for how the
            %   x-axis is formatted. Options are "linear", "log", and "dB".
            %   See QKDPlot for more details.
            % * yScaleStyle ("linear"): Name-value argument for how the
            %   y-axis is formatted. Options are "linear", "log", and "dB".
            %   See QKDPlot for more details.
            % * markerLineStyles ("x-"): Name-value argument for the marker
            %   and line style used.
            % * figAxis (axes(figure())): Axes object to plot on. Usefull
            %   when working with tiled/sub figures.
            %
            % See also QKDPlot
            arguments
                qkdInput (1,1) QKDSolverInput{QKDSolverInputMustHaveExactlyOneScanParam}
                results (:,1) struct
                options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle,["linear","log","dB"])} = "linear";
                options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle,["linear","log","dB"])} = "linear";
                options.markerLineStyles (1,1) string = "x-";
                options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure());
            end

            scanParamName = string(fieldnames(qkdInput.scanParameters));
            QKDPlot.plotParameters({results},scanParamName,...
                QKDPlot.keyRateTag,"xScaleStyle",options.xScaleStyle,...
                "yScaleStyle",options.yScaleStyle,...
                "markerLineStyles",options.markerLineStyles,...
                "figAxis",options.figAxis,"addLegend",false);   
        end

        function plotParametersFromFiles(fileNames,xParamName,yParamName,options)
            % plotParametersFromFiles Similar to QKDPlot.plotParameters but
            % uses a list of output files instead of a cell array of
            % results.
            %
            % Input:
            % * fileNames: String array of relative file path/names to plot
            %   results from.
            % * xParamName: Name of the parameter to plot on the x-axis.
            % * yParamName: Name of the parameter to plot on the y-axis.
            % * addLegend (true): Add a legend to the plot.
            % * legendNames (fileNames): Names used for each line sorted in
            %   the same order as the files. Defaults to the names of the
            %   files used.
            % * xScaleStyle ("linear"): Name-value argument for how the
            %   x-axis is formatted. Options are "linear", "log", and "dB".
            %   See QKDPlot for more details.
            % * yScaleStyle ("linear"): Name-value argument for how the
            %   y-axis is formatted. Options are "linear", "log", and "dB".
            %   See QKDPlot for more details.
            % * markerLineStyles (["x-","o-",".-","^-","+-",">-","s-","p-","*-"]): Name-value argument for the marker
            %   and line style used. If there are more lines than styles,
            %   then the system loops back to the start.
            % * figAxis (axes(figure())): Axes object to plot on. Usefull
            %   when working with tiled/sub figures.
            %
            % See also QKDPlot, QKDPlot.plotParameters
            arguments
                fileNames (:,1) string{mustBeFile}
                xParamName (1,1) string
                yParamName (1,1) string
                options.addLegend (1,1) logical =true;
                options.legendNames (:,1) string {mustBeEqualSize(options.legendNames,fileNames)} = fileNames;
                options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle,["linear","log","dB"])} = "linear";
                options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle,["linear","log","dB"])} = "linear";
                options.markerLineStyles (:,1) string = ["x-", "o-", ".-", "^-", "+-", ">-", "s-", "p-", "*-"];
                options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure());
            end

            % construct list of files to plot
            resultsSets = cell(size(fileNames));
            for index = 1:numel(fileNames)
                tempFile = load(fileNames(index));
                resultsSets{index} = tempFile.results;
            end

            QKDPlot.plotParameters(resultsSets,xParamName,yParamName,...
                "legendNames",options.legendNames,"xScaleStyle",options.xScaleStyle,...
                "yScaleStyle",options.yScaleStyle,"markerLineStyles",options.markerLineStyles,...
                "figAxis",options.figAxis, "addLegend",options.addLegend);
        end

        function plotParameters(resultsSets,xParamName,yParamName,options)
            % plotParameters plots the given parameters against each other
            % from the results struct. To support plotting multiple curves,
            % the results structs are stored together in a cell array.
            %
            % Input:
            % * resultsSets: Cell array where each element corresponds to a
            %   different run of the software to compare (Each element is a
            %   different struct array of results). For example to compare
            %   BB84 and SARG92 protocols, you would set this argument to
            %   {resultsBB84, resultsSARG92}.
            % * xParamName: Name of the parameter to plot on the x-axis.
            % * yParamName: Name of the parameter to plot on the y-axis.
            % * addLegend (true): Add a legend to the plot.
            % * legendNames (compose("data set %d",1:numel(resultsSets))): Names used for each line sorted in
            %   the same order as resultsSets. By default it uses a generic
            %   "data set x" name. 
            % * xScaleStyle ("linear"): Name-value argument for how the
            %   x-axis is formatted. Options are "linear", "log", and "dB".
            %   See QKDPlot for more details.
            % * yScaleStyle ("linear"): Name-value argument for how the
            %   y-axis is formatted. Options are "linear", "log", and "dB".
            %   See QKDPlot for more details.
            % * markerLineStyles (["x-","o-",".-","^-","+-",">-","s-","p-","*-"]): Name-value argument for the marker
            %   and line style used. If there are more lines than styles,
            %   then the system loops back to the start.
            % * figAxis (axes(figure())): Axes object to plot on. Usefull
            %   when working with tiled/sub figures.
            %
            % See also QKDPlot
            arguments
                resultsSets (:,1) cell {mustBeCellOf(resultsSets,"struct")}
                xParamName (1,1) string {mustBeANamedParameterOrKeyFlag(xParamName,resultsSets)}
                yParamName (1,1) string {mustBeANamedParameterOrKeyFlag(yParamName,resultsSets)}
                options.addLegend (1,1) logical =true;
                options.legendNames (:,1) string {mustBeEqualSize(options.legendNames,resultsSets)} = compose("data set %d",1:numel(resultsSets));
                options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle,["linear","log","dB"])} = "linear";
                options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle,["linear","log","dB"])} = "linear";
                options.markerLineStyles (:,1) string = ["x-", "o-", ".-", "^-", "+-", ">-", "s-", "p-", "*-"];
                options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure());
            end

            xUseDB = options.xScaleStyle == "dB";
            yUseDB = options.yScaleStyle == "dB";

            % extract and format data
            [xDataSets,simpleIndexingXFlag] = extractAndFormatDataSets(xParamName,resultsSets);
            [yDataSets,simpleIndexingYFlag] = extractAndFormatDataSets(yParamName,resultsSets);


            % dB scaling if requested

            % check to make sure we don't mix dB scaling with simple indexing
            if simpleIndexingXFlag && xUseDB || simpleIndexingYFlag && yUseDB
                throw(MException("plotParameters:CantSimpleIndexAndDBScale",...
                    "The data requires simple linear indexing which is not " + ...
                    "compatible with using dB scaling style."))
            end
            if xUseDB
                xDataSets = cellfun(@(x)-10*log10(x),xDataSets,"UniformOutput",false);
            end
            if yUseDB
                yDataSets = cellfun(@(x)-10*log10(x),yDataSets,"UniformOutput",false);
            end

            % label names
            xLabelName = formatLabelName(xParamName,simpleIndexingXFlag,xUseDB);
            yLabelName = formatLabelName(yParamName,simpleIndexingYFlag,yUseDB);

            % actually make the plot
            plotDataSets(options.figAxis,xDataSets,yDataSets,xLabelName,yLabelName,...
                options.xScaleStyle~="log",options.yScaleStyle~="log",...
                options.markerLineStyles,options.legendNames,options.addLegend);
        end
    end
end


%% extraction and formatting
function [dataSets, simpleIndexingFlag] = extractAndFormatDataSets(paramName,resultsSets)

dataSets = cell(size(resultsSets));
simpleIndexingFlag = false; %true means data isn't good for plotting use the index 1...n instead

for index =1:numel(resultsSets)
    dataSets{index} = extractParameter(paramName,resultsSets{index});

    %check if we need to simple index
    simpleIndexingFlag = ~all(cellfun(@(x)isscalar(x) && isreal(x), dataSets{index}));

    if simpleIndexingFlag
        break
    else
        dataSets{index} =cell2mat(dataSets{index});
    end
end

if simpleIndexingFlag
    % The simple indexing flag was triggered.
    % We can't use dB scaling and simple indexing at the same time

    % Make sure each set has the
    % same size and if so, convert to a simple indexing
    numelEach = cellfun(@numel,resultsSets);
    if all(numelEach == numelEach(1))
        dataSets{:} = 1:numelEach(1);
    else
        throw(MException("plotParameters:SimleAssignmentOnlyForSameSize",...
            "If the parameters can't be converted to the single real values," + ...
            " then each result must be the same size."))
    end
else

end
end

function labelName = formatLabelName(paramName,simpleIndexingFlag,useDB)

if paramName == QKDPlot.keyRateTag
    paramName = "key rate";
end

% handle _ character

paramName = replace(paramName,"_","\_");

% select name style
if simpleIndexingFlag
    labelName = paramName+" index";
elseif useDB
    labelName = paramName+" parameterized in dB";
else
    labelName = paramName;
end

end

function paramValues = extractParameter(paramName,results)
if paramName  == QKDPlot.keyRateTag
    paramValues = {results(:).keyRate};
else
    currentParams = [results(:).currentParams];
    paramValues = {currentParams.(paramName)};
end
end

%% main 1D plot function and helper tools
function plotDataSets(figAxis,xDataSets,yDataSets,xLabelName,yLabelName,linX,linY,markerLineStyles,legendNames, addLegend)
arguments
    figAxis (1,1) matlab.graphics.axis.Axes
    xDataSets (:,1) cell {mustBeCellOf(xDataSets,"double")}
    yDataSets (:,1) cell {mustBeCellOf(yDataSets,"double"),mustBeEqualSize(yDataSets,xDataSets)}
    xLabelName (1,1) string
    yLabelName (1,1) string
    linX (1,1) logical
    linY (1,1) logical
    markerLineStyles (:,1) string
    legendNames (:,1) string
    addLegend (1,1) logical
end

if linX
    figAxis.XScale = 'linear';
else
    figAxis.XScale = 'log';
end
if linY
    figAxis.YScale = 'linear';
else
    figAxis.YScale = 'log';
end

figAxis.XLabel.String = xLabelName;
figAxis.YLabel.String = yLabelName;

for index = 1:numel(xDataSets)
    if index == 1
        hold(figAxis,"on");
    end
    markerLineStyle = markerLineStyles(mod(index-1,numel(markerLineStyles))+1);
    plot(figAxis,xDataSets{index},yDataSets{index},markerLineStyle,...
        "DisplayName",legendNames(index));
end
hold(figAxis,"off");

if addLegend
    legend(figAxis)
end

end

%% validation functions
function mustBeANamedParameterOrKeyFlag(paramName,resultsSets)
if paramName == QKDPlot.keyRateTag
    return
end
if ~all(cellfun(@(x)ismember(paramName,string(fieldnames(x(1).currentParams))), resultsSets))
    throwAsCaller(MException("plotParameters:NotAParameter",...
        sprintf("%s was not a parameter name used in all the results sets " + ...
        "or the key rate tag,'%s'.",paramName,QKDPlot.keyRateTag)))
end
end

function QKDSolverInputMustHaveExactlyOneScanParam(qkdInput)
numParams = cellfun(@numel,{qkdInput.scanParameters});
if any(numParams~=1)
    throw(MException("plotResults:MustHaveExactlyOneScanParam",...
        "The QKDSolverInput must have exaclty one scan parameter for plotting."))
end
end