classdef PatternSearchAlgorithm
    properties
        pollString
        isNUPS
    end

    methods
        function alg = PatternSearchAlgorithm(pollString, isNUPS)
            alg.pollString = pollString;
            alg.isNUPS = isNUPS;
        end
    end

    enumeration
        GPS2N ("GPSPositiveBasis2N", 0)
        GPSNp1 ("GPSPositiveBasisNp1", 0)
        GSS2N ("GSSPositiveBasis2N", 0)
        GSSNp1 ("GSSPositiveBasisNp1", 0)
        MADS2N ("MADSPositiveBasis2N", 0)
        MADSNp1 ("MADSPositiveBasisNp1", 0)
        OrthoMADS2N ("OrthoMADSPositiveBasis2N", 0)
        OrthoMADSNp1 ("OrthoMADSPositiveBasisNp1", 0)
        NUPS ("nups", 1)
        NUPSGPS ("nups-gps", 1)
        NUPSMADS ("nups-mads", 1)
    end
end
