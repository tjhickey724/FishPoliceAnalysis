function [ ] = gandsfigs(e, cmd )
%gandsfigs generates the figures for the Games and Simulation paper
%   Detailed explanation goes here

%auditory = [4,5,6,8,11,13,14,21,22,23,24,26,27,30,33,34,36,39,40,41];
%visual = [1,2,3,7,9,10,12,15,16,17,18,19,20,25,28,29,31,32,35,37,38,42,43];

auditory = [4,5,6,8,11,13,14,21,22,23,24,26,27,30,33,34,36]; %novices
visual = [1,2,3,7,9,10,12,15,16,17,18,19,20,25,31,32,35,37]; %novices
pros = [39,40,41,28,29,38,42,43];
doubles = [26,41,25,43,17,38,27,40]
triples = [9 28 39 11 29 42];
    switch cmd
        case 'fig7x7'
            figure(1)
            subplot(1,1,1)
            for i=[1:length(auditory)]
                subplot(5,4,i);
                runDemo(e,'gameplay4a1',auditory(i))
            end
            figure(2)
            subplot(1,1,1)
            for i=[1:length(visual)]
                subplot(5,4,i);
                runDemo(e,'gameplay4a1',visual(i))
            end
            figure(3)
            subplot(1,1,1)
            for i=[1:length(pros)]
                subplot(5,4,i);
                runDemo(e,'gameplay4a1',pros(i));
            end
            figure(4)
            subplot(1,1,1)
            for i = [1:6]
                subplot(4,3,i);
                runDemo(e,'gameplay4a1',triples(i));
            end
            figure(5)
            for i = [1:8]
                subplot(4,2,i);
                runDemo(e,'gameplay4a1',doubles(i));
            end

        case 'e1plot'
            for i=[1:15]
                subplot(4,4,i);
                runDemo(e,'gameplay4',i);
                title(i)
            end
        case 'e2plot'
            for i=[16:29]
                subplot(4,4,i-15);
                runDemo(e,'gameplay4',i);
                title(i)
            end
        case 'e3plot'
            for i=[30:43]
                subplot(4,4,i-29);
                runDemo(e,'gameplay4',i);
                title(i)
            end
            
    end
end

