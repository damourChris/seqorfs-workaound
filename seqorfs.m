function [orfOutput,seqLen] = seqorfs(seq,varargin)
% 提取序列的-3~+3orf框的位置信息
%SEQSHOWORFS highlights the open reading frames (ORFs) in a sequence.
%
%   SEQSHOWORFS(SEQ) displays the sequence SEQ in the MATLAB figure window,
%   with all ORFs highlighted. SEQSHOWORFS returns a structure of the start
%   and stop positions of the ORFs in each reading frame. The Standard
%   genetic code is used with start codon 'AUG' and stop codons
%   'UAA','UAG','UGA'.
%
%   SEQSHOWORFS(...,'FRAMES',READINGFRAMES) specifies the reading frames to
%   display. READINGFRAMES can be any of 1,2,3,-1,-2,-3. Frames -1, -2, and
%   -3 correspond to the first, second, and third reading frames of the
%   reverse complement of SEQ. To display multiple frames at once, use a
%   vector of the frames to be displayed, or use 'all' to show all frames
%   at once. The default is [1 2 3].
%
%   SEQSHOWORFS(...,'GENETICCODE',CODE) specifies the genetic code to be
%   used for finding open reading frames. CODE can be an integer specifying
%   a code ID, a character vector or string specifying a code name, or a
%   structure created using the function GENETICCODE. See help for
%   geneticcode for the full list of supported IDs and names.
%
%   SEQSHOWORFS(...,'MINIMUMLENGTH',L) sets the minimum number of codons
%   for an ORF to be considered valid. The default is 10.
%
%   SEQSHOWORFS(...,'ALTERNATIVESTARTCODONS',true) uses alternative start
%   codons. For example in the Human Mitochondrial Genetic Code, AUA and
%   AUU are known to be alternative start codons.
%
%   See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
%   for more details of alternative start codons.
%
%   SEQSHOWORFS(...,'COLOR',COLOR) selects the color used to highlight the
%   open reading frames in the output display. The default color scheme is
%   blue for the first reading frame, red for the second, and green for the
%   third frame. COLOR can be a 1x3 RGB vector whose elements specify the
%   intensities of the red, green and blue component of the color; the
%   intensities can be in the range [0 1]. COLOR can also be a character
%   from the following list:
%
%            'b'     Blue
%            'g'     Green
%            'r'     Red
%            'c'     Cyan
%            'm'     Magenta
%            'y'     Yellow
%
%   To specify different colors for the three reading frames, use a 1x3
%   cell array of color values. If you are displaying reverse complement
%   reading frames, then COLOR should be a 1x6 cell array of color values.
%
%   SEQSHOWORFS(...,'COLUMNS',COLS) specifies how many columns per line to
%   use in the output. The default is 64.
%
%   Example:
%
%       HLA_DQB1 = getgenbank('NM_002123');
%       seqshoworfs(HLA_DQB1.Sequence);
%
%   See also CODONCOUNT, CPGISLAND, GENETICCODE, REGEXP, SEQDISP,
%   SEQSHOWWORDS, SEQVIEWER, SEQWORDCOUNT.

%   SEQSHOWORFS(...,'TRANSLATE',TF) adds translation for each ORF to the output
%   structure array when set TRUE. The default is FALSE.

%   Copyright 2002-2012 The MathWorks, Inc.


% set up some defaults

if nargin > 0
    seq = convertStringsToChars(seq);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

wrap = 64;
theFrames = 1:3;
revComp = false;
startCodons = {'ATG','CTG','TTG'}; % Defaults to standard genetic code
stopCodons = {'TAA','TAG','TGA'}; % Defaults to standard genetic code
minLength = 10;

colors = {'0000FF','FF0000','00BB00','0000FF','FF0000','00BB00'};
alternativestartcodons = false;
noDisplay = true;
transFlag = false;
gcode = 1;

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

% deal with the various inputs
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:seqshoworfs:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'color','columns','geneticcode','alternativestartcodons',...
        'frames','nodisplay','minimumlength','translate'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:seqshoworfs:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqshoworfs:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 %color
                    if ~iscell(pval)
                        pval = {pval,pval,pval,pval,pval,pval}; %#ok<AGROW>
                    end

                    for count = 1:length(pval)
                        val = pval{count};
                        colors{count} = setcolorpref(val, colors{count});
                    end
                case 2% wrap
                    wrap = pval;
                case 3 % genetic code
                    gcode = pval;
                    if isstruct(pval)
                        numFields = numel(fieldnames(pval));
                        % allow revgeneticcode created struct
                        if numFields == 23  && isfield(pval, 'Name')
                            pval = geneticcode(pval.Name);
                            numFields = numel(fieldnames(pval));
                        end
                        if numFields ~= 66 || ~isfield(pval, 'Name')
                            error(message('bioinfo:seqshoworfs:BadCodeStruct'));
                        end
                        gc = pval;
                    else
                        try
                            gc = geneticcode(pval);
                        catch allExceptions
                            if isnumeric(pval)
                                error(message('bioinfo:seqshoworfs:BadGeneticCodeID', num2str( pval )))
                            else
                                error(message('bioinfo:seqshoworfs:BadGeneticCode', pval))
                            end
                        end
                    end
                    startCodons = gc.Starts;
                    fn = fieldnames(gc);
                    stopCodons = fn(strcmp(struct2cell(gc),'*'))';

                case 4 % use alternative start codons
                    alternativestartcodons = bioinfoprivate.opttf(pval);
                    if isempty(alternativestartcodons)
                        error(message('bioinfo:seqshoworfs:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 5  %frames
                    if isnumeric(pval)
                        theFrames = double(pval);
                        if ~all(ismember(theFrames,[1,2,3,-1,-2,-3]))
                            error(message('bioinfo:seqshoworfs:InvalidReadingFrameNum'));
                        end
                    elseif ischar(pval)
                        if strcmpi(pval,'all')
                            theFrames = [1 2 3 -1 -2 -3];
                        else
                            error(message('bioinfo:seqshoworfs:InvalidReadingFrameChar'));
                        end
                    end
                    if any(theFrames < 0) % need rcomplement
                        revComp = true;
                        theFrames(theFrames<0) = 3-theFrames(theFrames<0);
                    end
                case 6 % nodisplay -- use this for testing with second output
                    noDisplay = bioinfoprivate.opttf(pval);
                    if isempty(noDisplay)
                        error(message('bioinfo:seqshoworfs:DisplayInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 7 % minimum length
                    if isnumeric(pval) && pval > 0
                        minLength = pval;
                    else
                        error(message('bioinfo:seqshoworfs:MinimumLengthNotPositiveInt'));
                    end
                case 8 % translate flag -- use this for seqviewer
                    transFlag = bioinfoprivate.opttf(pval);
                    if isempty(transFlag)
                        error(message('bioinfo:seqshoworfs:TranslateOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

seqLen = length(seq);

if ~noDisplay && (seqLen * length(theFrames) > 1000000)  % large display will run out of Java Heap Space
    error(message('bioinfo:seqshoworfs:ShowORFSLimit'))
end

% set up output
orfOutput.Start = [];
orfOutput.Stop = [];
if transFlag
    orfOutput.Translation = {''};
end
orfOutput = repmat(orfOutput,1,max(theFrames));

% set default startCodons
if alternativestartcodons == false
    startCodons = {'ATG'};
end

% work also with RNA sequences
stopCodons = strrep(stopCodons,'T','[TU]');
startCodons = strrep(startCodons,'T','[TU]');

% if we need to do reverse complement then make a copy
if revComp
    rseq = seqrcomplement(seq);
    rstarts = sort(cell2mat(regexpi(rseq,startCodons)));
    rstops = sort(cell2mat(regexpi(rseq,stopCodons)));
    fstarts = sort(cell2mat(regexpi(seq,startCodons)));
    fstops = sort(cell2mat(regexpi(seq,stopCodons)));
    fseq = seq;
else
    starts = sort(cell2mat(regexpi(seq,startCodons)));
    stops = sort(cell2mat(regexpi(seq,stopCodons)));
    clear fseq;
end

revSeq = [];
for frame = theFrames
    if revComp  % need to deal with switching between two directions
        if frame < 4
            starts = fstarts;
            stops = fstops;
            revSeq = fseq;
        else
            starts = rstarts;
            stops = rstops;
            revSeq = rseq;
        end
    end

    % We now have to deal with figuring out which starts and stops are in
    % this frame.
    framestartcodons = starts(rem(starts,3)==rem(frame,3));
    framestopcodons = stops(rem(stops,3)==rem(frame,3));

    % and we only care about stops after a start and starts after a stop
    % so we need to throw away any stops not in ORF and any starts inside
    % an already open reading frame.
    stopsAndStart = union(framestartcodons,framestopcodons);
    [~,StartIndices ] = ismember(framestartcodons,stopsAndStart);
    [~,StopIndices ] = ismember(framestopcodons,stopsAndStart);
    StartDiffs = diff([-1 StartIndices]);
    StopDiffs = diff([0 StopIndices]);

    frameStart = framestartcodons(StartDiffs ~= 1);
    frameStop = framestopcodons(StopDiffs ~= 1);

    % look for small ORFs and remove them
    if ~isempty(frameStop)
        smallORFs = find((frameStop-frameStart(1:numel(frameStop)))< minLength*3);
        frameStart(smallORFs) = [];
        frameStop(smallORFs) = [];
    end
    % clean up empty results
    if isempty(frameStart)
        frameStart = [];
    end
    if isempty(frameStop)
        frameStop = [];
    end

    % add on an extra holding place just to make sure that frameStop is as
    % long as frameStart.
    numOrf = numel(frameStart);
    if numOrf > numel(frameStop)
        frameStop(end+1) = seqLen+1; %#ok
    end

    % save the output
    orfOutput(frame).Start = frameStart;
    orfOutput(frame).Stop = frameStop;

    if isempty(frameStart)
        frameStart=seqLen + 1;
    end
end


orfsSegs = struct('length',[],'shift',[],'iscomplete',[],'x',[],'y',[]);
k = 1;

num_frames = length(orfOutput);

for i = 1:num_frames
    ose = orfOutput(i);
    segLens = ose.Stop-ose.Start;
    j = 1;
    shiftc = i;
    for Segrec = segLens
        orfsSegs(k).x = ose.Start(j);
        orfsSegs(k).y = ose.Stop(j);
        if i>3
            shiftc = 3 - i;
            orfsSegs(k).y = seqLen-ose.Stop(j);
            orfsSegs(k).x = seqLen-ose.Start(j);
        end
        orfsSegs(k).shift =  shiftc;
        if ose.Stop(j) > seqLen
            orfsSegs(k).iscomplete = 0;
        else
            orfsSegs(k).iscomplete = 1;
        end
        orfsSegs(k).length = Segrec;

        j = j + 1;
        k = k + 1;
    end

end % end for

orfOutput = struct2table(orfsSegs);
orfOutput.Properties.VariableNames = {'length','shift','iscomplete','start','stop'};
end % function
