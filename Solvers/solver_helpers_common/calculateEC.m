%% Function name: calculateEC
% This code calculates H(X|Y), H(Y|X) and IXY from a given probability
% distribution of outcomes. 
% 
% Please input a normalized probability distribution table for the
% calculation.

function [HX_Y, HY_X, IXY] = calculateEC(prob_dist)   	
   
   
    [nRows,nColumns] = size(prob_dist);
   
    px_dist = sum(prob_dist,2);
    py_dist = sum(prob_dist);
    HXY = 0;
    HX = 0;
	for iRow =1:nRows
        for jColumn=1:nColumns
            if prob_dist(iRow,jColumn)~=0 
                HXY = HXY - prob_dist(iRow,jColumn) * log2(prob_dist(iRow,jColumn));
            end
            
        end
        if px_dist(iRow) ~=0
           HX = HX - px_dist(iRow) * log2(px_dist(iRow)); 
        end
        
    end
 	HY = 0;
	for jColumn =1:nColumns
        if py_dist(jColumn) ~=0
            HY = HY - py_dist(jColumn) * log2(py_dist(jColumn)); 
        end
	end
    HY_X = HXY - HX;
    HX_Y = HXY - HY;
    IXY  = HY - HY_X;
end