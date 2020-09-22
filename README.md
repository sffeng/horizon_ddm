# horizon_ddm
Code to reproduce figures and analysis for Feng, Wang, Zarnescu, and Wilson: "The dynamics of explore-exploit decisions reveal a signal-to-noise mechanism for random exploration"


TODO: Add comment on the notational change in the way the code is written vs the way the paper is written.

In the MATLAB code, the notation for the DDM is slightly different from that used in the paper, following instead that Bogacz et al 2006.  The differences are in the threshold and bias (starting point).  In the paper, the starting point was \alpha \beta, whereas in Bogacz et al this is treated as a single separate parameter, X0. In the paper, the upper/lower thresholds are at \beta and 0 (implying that 0 <= \alpha \beta <= \beta), whereas in Bogacz et al, thresholds are between z and -z for some positive parameter z. 


Bogacz, R., Brown, E., Moehlis, J., Holmes, P., & Cohen, J. D. (2006). The physics of optimal decision making: a formal analysis of models of performance in two-alternative forced-choice tasks. Psychological review, 113(4), 700.
