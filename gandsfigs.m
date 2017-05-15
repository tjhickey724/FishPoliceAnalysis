function [ ] = gandsfigs(e, cmd )
%gandsfigs generates the figures for the Games and Simulation paper
%   Detailed explanation goes here

auditory = [4,5,6,8,11,13,14,21,22,23,24,26,27,30,33,34,36,39,40,41];
visual = [1,2,3,7,9,10,12,15,16,17,18,19,20,25,28,29,31,32,35,37,38,42,43];

    switch cmd
        case 'fig7x7'
            figure(1)
            subplot(1,1,1)
            for i=[1:length(auditory)]
                subplot(5,5,i);
                runDemo(e,'gameplay4a1',auditory(i))
            end
            figure(2)
            subplot(1,1,1)
            for i=[1:length(visual)]
                subplot(5,5,i);
                runDemo(e,'gameplay4a1',visual(i))
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

