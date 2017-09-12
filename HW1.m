% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
%x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
 x = 3; y = '5'; %mixed

%your code goes here
if isnumeric(x) == false
    x = str2num(x);
end
if isnumeric(y) == false
    y = str2num(y);
end
%output your answer
answer = x+y;
disp(answer)
disp ('Terminated Analysis')%Aids in reading Command window runs.

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with 
% N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 
clear all
N = 500; % define sequence length

randNum_seq = randi(4,[1,N]) %generates an array of numbers 1-4.
%this will correspond to each of the 4 nucleotide bases.

ii = 0;
while ii < N
    ii = ii+1;% counter starts/increases
    switch randNum_seq(ii)
        case 1
            rand_seq(1,ii) = 'A';
        case 2
            rand_seq(1,ii )= 'T';
        case 3
            rand_seq(1,ii) = 'C';
        case 4
            rand_seq(1,ii) = 'G';
    end

end
disp(rand_seq)
disp ('Terminated Analysis')%Aids in reading Command window runs.


%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

StartCodon = ['ATG'];% Defines starting codon text
StopCodon1 = ['TAA'];% Defines Stopping codon text
StopCodon2 = ['TGA'];% Defines Stopping codon text
StopCodon3 = ['TAG'];% Defines Stopping codon text

Start = strfind(rand_seq,StartCodon); %returns an array of all Start codons
Stop1 = strfind(rand_seq,StopCodon1); %returns an array of all stop codons for TAA
Stop2 = strfind(rand_seq,StopCodon2); %returns an array of all stop codons for TGA
Stop3 = strfind(rand_seq,StopCodon3); %returns an array of all stop codons for TAG
Stop = [Stop1,Stop2,Stop3]; Stop = sort(Stop)%stop codon listed in numerical order 

q = 1;
DistMax = -1;

while q<=length(Start)
    checkdist = Stop-Start(1,q); %checks distance between every stop and 1 element of start
    p = 1;
    for ii = 1:length(checkdist) %finds the first, therefore shortest, 
                                 %positive distance. 
                                 %All subsequent distances will be
                                 %terminated by this first distance.
        if sign(checkdist(1,p)) == 1
            ArrayMax = checkdist(1,p)
            break
            break
        else
            p = p+1;
        end
    end
        
    if ArrayMax > DistMax %Saves only the longest distance. 
        DistMax = ArrayMax 
    end
    q = q+1 %Counter
end
DistMax
disp ('Terminated Analysis') %Aids in reading Command window runs.

%%
%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.
clear all
RunNumber = 1000;
probmat = 0;

for z = 1:RunNumber
    
N = 500; % define sequence length
randNum_seq = randi(4,[1,N]); %generates an array of numbers 1-4.
%this will correspond to each of the 4 nucleotide bases.

ii = 0;
while ii < N
    ii = ii+1;% counter starts/increases
    switch randNum_seq(ii)
        case 1
            rand_seq(1,ii) = 'A';
        case 2
            rand_seq(1,ii )= 'T';
        case 3
            rand_seq(1,ii) = 'C';
        case 4
            rand_seq(1,ii) = 'G';
    end
end %Generates the random DNA sequence

StartCodon = ['ATG'];% Defines starting codon text
StopCodon1 = ['TAA'];% Defines Stopping codon text
StopCodon2 = ['TGA'];% Defines Stopping codon text
StopCodon3 = ['TAG'];% Defines Stopping codon text
Start = strfind(rand_seq,StartCodon); %returns an array of all Start codons
Stop1 = strfind(rand_seq,StopCodon1); %returns an array of all stop codons for TAA
Stop2 = strfind(rand_seq,StopCodon2); %returns an array of all stop codons for TGA
Stop3 = strfind(rand_seq,StopCodon3); %returns an array of all stop codons for TAG
Stop = [Stop1,Stop2,Stop3]; Stop = sort(Stop);%stop codon listed in numerical order 

q = 1;
DistMax = -1;
ArrayMax = 0;
while q<=length(Start)
    checkdist = Stop-Start(1,q); %checks distance between every stop and 1 element of start
    p = 1;
    for ii = 1:length(checkdist) %finds the first, therefore shortest, 
                                 %positive distance. 
                                 %All subsequent distances will be
                                 %terminated by this first distance.
        if sign(checkdist(1,p)) == 1
            ArrayMax = checkdist(1,p);
            break
            break
        else
            p = p+1;
        end
    end
        
    if ArrayMax > DistMax %Saves only the longest distance. 
        DistMax = ArrayMax ;
    end
    q = q+1; %Counter
end %Finds the max ORF distance

if DistMax >= 50 %finds the number of distances > than 50bp
    probmat = probmat+1;
end

end
probmat/RunNumber;

disp('the probability is:')
disp(probmat/RunNumber)

disp ('Terminated Analysis') %Aids in reading Command window runs.
%%
%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length.
clear all
RunNumber = 1000;
probmat = 0;
N = 58; % define sequence length

w = 1;
for iii = 1:N
    
for z = 1:RunNumber
randNum_seq = randi(4,[1,N]); %generates an array of numbers 1-4.
%this will correspond to each of the 4 nucleotide bases.

ii = 0;
while ii < N
    ii = ii+1;% counter starts/increases
    switch randNum_seq(ii)
        case 1
            rand_seq(1,ii) = 'A';
        case 2
            rand_seq(1,ii )= 'T';
        case 3
            rand_seq(1,ii) = 'C';
        case 4
            rand_seq(1,ii) = 'G';
    end
end %Generates the random DNA sequence

StartCodon = ['ATG'];% Defines starting codon text
StopCodon1 = ['TAA'];% Defines Stopping codon text
StopCodon2 = ['TGA'];% Defines Stopping codon text
StopCodon3 = ['TAG'];% Defines Stopping codon text
Start = strfind(rand_seq,StartCodon); %returns an array of all Start codons
Stop1 = strfind(rand_seq,StopCodon1); %returns an array of all stop codons for TAA
Stop2 = strfind(rand_seq,StopCodon2); %returns an array of all stop codons for TGA
Stop3 = strfind(rand_seq,StopCodon3); %returns an array of all stop codons for TAG
Stop = [Stop1,Stop2,Stop3]; Stop = sort(Stop);%stop codon listed in numerical order 

q = 1;
DistMax = -1;
ArrayMax = 0;
while q<=length(Start)
    checkdist = Stop-Start(1,q); %checks distance between every stop and 1 element of start
    p = 1;
    for ii = 1:length(checkdist) %finds the first, therefore shortest, 
                                 %positive distance. 
                                 %All subsequent distances will be
                                 %terminated by this first distance.
        if sign(checkdist(1,p)) == 1
            ArrayMax = checkdist(1,p);
            break
            break
        else
            p = p+1;
        end
    end
        
    if ArrayMax > DistMax %Saves only the longest distance. 
        DistMax = ArrayMax ;
    end
    q = q+1; %Counter
end %Finds the max ORF distance

if DistMax >= 50 %finds the number of distances > than 50bp
    probmat = probmat+1;
end

end %Determins the probablity of an ORF > 50bp.
Probability (w,:) = probmat/RunNumber;

w=w+1

end
plot(Probability)
disp ('Terminated Analysis') %Aids in reading Command window runs.

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

%I believe that there should be a linear relationship between the number of
%ORFs > 50bp and the size of the sequence. This is reflected in my plot.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 
clear all

fid = fopen('qPCRdata.txt', 'r');

countline = 3;
for i = 1:countline
line1 = fgetl(fid);
end

q = 1;
stopline = 'G1';
while line1 ~= -1
    var1 = strsplit(line1);
    str_mat(q,:) = var1;
    if strcmp(str_mat(q,5),'72')
        break
    end
    line1 = fgetl(fid);
    q = q+1;
end

num_mat = [str_mat(:,3),str_mat(:,6)]


% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc.

num_mat = num_mat';
q = 1;
p = 1;
for ii = 1:6 %makes rows 1-6 of 12 columns
    wellplate(ii,:) = num_mat(2,q:q+11)
    p = p+1;
    q = q+12;
end
size(wellplate)%checks array size
wellplate = cellfun(@str2num,wellplate)%changes array data type to num.
disp('Terminated Analysis')

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

qq = 1
pp = 1
for ix = 1:6
    y(pp,:) = 2.^(wellplate (1,1:9)- wellplate(qq,1:9))
    diff_mat(pp,:) = wellplate (1,10:12)-wellplate(qq,10:12)
    qq = qq+1;
    pp = pp+1;
end

nextcol = 1;

    for iiix = 1:3
    final_mat(:,nextcol) = y(:,nextcol)- diff_mat(:,nextcol)
    nextcol = nextcol+3;
    end



%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty