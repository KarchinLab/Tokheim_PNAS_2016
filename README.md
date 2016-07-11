# Cancer driver evaluation

A major goal of the huge public investment in large-scale cancer sequencing has been to find the majority of driver genes. Robust computational prediction of drivers from small somatic variants is critical to this mission, and it is essential that the best methods be identified. While many such methods have been proposed, it has been difficult to evaluate them because there is no gold standard to use as a benchmark. Here we developed an evaluation framework for driver gene prediction methods that does not require a gold standard. The framework includes a large set of small somatic mutations from a wide range of cancer types and five evaluation metrics. We propose it can be used to systematically evaluate new prediction methods and compare them to existing methods. The evaluations include:

* Fraction overlap with Cancer Gene Census (CGC)
* Method concensus
* P-value distribution
* Number of predicted driver genes
* Prediction consistency

## Jupyter notebooks

We have prepared our evaluation framework into jupyter notebooks (.ipynb files). You can either view
them on github or execute the evaluation if you install [jupyter](http://jupyter.org/).

## Installation

You will need the following python packages if you wish to execute the code:

* numpy
* scipy
* matplotlib
* seaborn
* pandas
