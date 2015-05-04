# GrangerCausality
This library is the java implementation of C.W.J. Grangers popular causality indicator metric, also called [Granger-Causality][1]. The implementation follows his paper from 1969 \[1\]. An easy description of the Granger-Causality would be as follows: "Does the historical information of varibale *x* help in predicting variable *y*?"  
There are other open source java implementations of the Granger-Causality, e.g., Sergey Edunov's GNU implementation [3]. However, these implementations are focussed on a bivariate space, restricting the user to 2 variables (x and y) only. In most cases, the universe that is being analysed consists of more varibles (is multivariate). The more appropriate question in such cases is therefore: "Does the historical information of varibale *x* help in predicting variable *y* in a given universe *u*?" This requires the extension of the bivariate Granger-Causality to a multivariate space, which is provided on top of the bivariate implementation with this library.  

[1]: https://en.wikipedia.org/wiki/Granger_causality  
[3]: https://code.google.com/p/jquant/source/browse/trunk/src/ru/algorithmist/jquant/math/GrangerTest.java?r=8  

\[1\] Granger, C. W. J. (1969). *"Investigating Causal Relations by Econometric Models and Cross-spectral Methods"*. Econometrica 37 (3): 424–438.  

