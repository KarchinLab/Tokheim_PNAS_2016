# Cancer driver evaluation

A major goal of the huge public investment in large-scale cancer sequencing has been to find the majority of driver genes. Robust computational prediction of drivers from small somatic variants is critical to this mission, and it is essential that the best methods be identified. While many such methods have been proposed, it has been difficult to evaluate them because there is no gold standard to use as a benchmark. Here we developed an evaluation framework for driver gene prediction methods that does not require a gold standard. The framework includes a large set of small somatic mutations from a wide range of cancer types and five evaluation metrics. We propose it can be used to systematically evaluate new prediction methods and compare them to existing methods. The evaluations include:

1. Fraction overlap of predicted driver genes with Cancer Gene Census (CGC)
2. Method concensus
3. P-value distribution
4. Prediction consistency
5. Number of predicted driver genes

The number of predicted driver genes can be informative along side other factors to understand whether there is an inflation of false positive driver genes or particularly low power. Due to both extremes being unfavorable, we do not include it as a measurement for ranking each method. The average ranking of a method on measures 1-4 can be used to assess relative performance.

## Jupyter notebooks

We have prepared our evaluation framework into jupyter notebooks (.ipynb files). You can either view
them on github or execute the evaluation if you install [jupyter](http://jupyter.org/). We recommend that you first
view the [Introduction.ipynb](https://github.com/KarchinLab/CancerDriverGeneEvaluation/blob/master/Introduction.ipynb) file for more details.

## Installation

You will need the following python packages if you wish to execute the python code in Jupyter:

* numpy
* scipy
* matplotlib
* seaborn
* pandas

## Citation

Collin J. Tokheim, Nickolas Papadopoulos, Kenneth W. Kinzler, Bert Vogelstein, and Rachel Karchin. Evaluating the evaluation of cancer driver genes. PNAS 2016 ; published ahead of print November 22, 2016, doi:10.1073/pnas.1616440113
