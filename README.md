# DAG-BetaScores

The function now distinguishes between two cases for each potential parent node:
- If the parent node is already present in the DAG, it calculates the beta score as
  the difference between the score with that parent and other existing parents, and the score without that specific parent.
  
- If the parent node is not present in the DAG, it calculates the beta score as
  the difference between the score with the parent node added and the score without the parent node.
