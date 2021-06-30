# TSNet

TSNet (Tumor Specific Net) is a new method which constructs tumor-cell specific gene/protein co-expression networks based on gene/protein expression profiles of tumor tissues. TSNet treats the observed expression profile as a mixture of expressions from different cell types and explicitly models tumor purity percentage in each tumor sample.

Usage
The main function - deNet - takes as input gene or protein expression profiles of bulk tumor data. This function returns the estimated purity, mean expression of different markers in tumor cells and non-tumor cells (combination of immune and stromal cells) and co-expression networks in tumor cells and non-tumor cells.
