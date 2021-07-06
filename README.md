# codmi
<b>codmi</b> is an R algorithm which allows for the inclusion of oncological COVID-19 deaths in traditional survival analysis.<br>
For a given sample of time-to-event data, typically oncological deaths with censored observations, where covid deaths are also observed, <b>codmi</b> algorithm includes covid-death cases in the standard data by mean-imputation. Data completion is done through an expectation-maximization algorithm based on the Kaplan-Meier estimator. The completed data are equipped with the corresponding Kaplan-Meier survival function estimates and the appropriate confidence intervals. These are based on an extended Greenwood’s formula where the variance is adjusted for the mean-imputation.

<br>
This branch of the repository contains the basic functions of the algorithm and a suite of examples.
A self-contained installable package is in the works.

# usage
* <b>codmi</b> needs a working installation of the [R programming language](http://www.r-project.org/) to run;
* install the ```survival``` and the ```stats``` third party packages;
* download the repository and execute the code in ```examples.R```
