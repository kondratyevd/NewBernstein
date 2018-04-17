# Problem with RooBernstein()
In RooFit, fit with Bernstein polynomials (RooBernstein() class) doesn't work correctly in multiple subranges.

The problem occurs during evaluation of the integral of a polynomial. Bernstein polynomials are defined in range [0, 1]. The original algorithm, contracts any range down to [0, 1], then evaluates a definite integral of the polynomial over [0, 1] and then rescales the function back. 

That obviously leads to incorrect behaviour when dealing with multiple subranges.

# Solution
**NewBernstein()** is designed using the approach similar to that described here: https://sft.its.cern.ch/jira/browse/ROOT-6664 .

The following is implemented: the range [0, 1] is mapped to all of the input subranges in a way that the minimal *x* value among all subranges is mapped to 0 and the maximum value is mapped to 1. All the other ends of subranges end up somwhere between 0 and 1. 

For example, two subranges [100, 125] and [175, 200] will be mapped to [0, 0.25] and [0.75, 1].

After that, the integrals over all the subranges are evaluated and the function is rescaled back.


# Implementation in ROOT:

Follow the directions here ("Extending ROOT with custom classes"): 
http://webhome.phy.duke.edu/~dmb60/the-guide/
