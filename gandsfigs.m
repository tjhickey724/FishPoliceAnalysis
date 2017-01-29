function [ ] = gandsfigs(e, cmd )
%gandsfigs generates the figures for the Games and Simulation paper
%   Detailed explanation goes here


    switch cmd
        case 'fig7x7'
            for i=[1:43]
                subplot(7,7,i);
                runDemo(e,'gameplay4',i)
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

