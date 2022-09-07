# Brain segmentation

[Presentation](https://docs.google.com/presentation/d/1RexJKJhKoz3l7kB0DZQ_h-yEBJ2Py50UVGyMlsiq0cQ/edit?usp=sharing), but it is on russian

| Algorithm      |  MEL1  |  MEL2  |  MEL3  |  SIM2  |  SIM3  |  SIM4  |
| :---           | :----: | :----: | :----: | :----: | :----: | :----: |
| Best classical | 0.922  | 0.925  | 0.909  | 0.925  | 0.949  | 0.926  |
| Genetic        | 0.906  | 0.926  | 0.927  | 0.918  | 0.892  | 0.94   |

* MEL and SIM are different types of flies
* Best classical is one of these algorithm(with best score for each brain): [K-means](https://scikit-learn.org/stable/modules/clustering.html#k-means), [Agglomerative clustering](https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering), [BIRCH](https://scikit-learn.org/stable/modules/clustering.html#birch)
 
