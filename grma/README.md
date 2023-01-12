# grma
> grma - Graph based Matching

grma is a package for finding HLA matches using graphs approach.
The matching is based on [grim](https://github.com/nmdp-bioinformatics/py-graph-imputation)'s imputation.

## Building The Donors' Graph

The donors' graph is a graph implemented using a LOL (List of Lists) representation written in cython for better
time and memory efficiency.

Before building the donors' graph, all the donors' HLAs must be imputed using `grim`.
Then all the imputation files must be saved under the same directory.

```python
from grma.donorsgraph.build_donors_graph import BuildMatchingGraph

build_matching = BuildMatchingGraph("path/to/donors/directory")
graph = build_matching.graph  # get the donors' graph
build_matching.to_pickle("path/to/save.pkl")  # save the donors' graph to pickle
```

## Search & Match

Find matches up to 3 mismatches and get a `pandas.DataFrame` object of the matches sorted by number of mismatches and score.

### Set Database
In order to get in the matching results more information about the donors than the matching information,
one can set a database that has all the donors' information in it.
The database must be a `pandas.DataFrame` that its indexes are the donors' IDs.

```python
import pandas as pd
from grma.match import set_Database

donors = [0, 1, 2]
database = pd.DataFrame([[30], [32], [25]], columns=["Age"], index=donors)

set_Database(database)
```

Note: `set_Database()` need to be called only once before
performing the matching if additional fields are ment to be added to the results.

### Matching

Get the donors' graph from saved pickle file:
```python
from grma.match import Graph
graph = Graph.from_pickle("path/to/save.pkl")
```

The matching result is a dictionary that maps patients IDs to a `pandas.DataFrame` of the matched donors.
It gets:
 - A path to the file of the patients' typing.
 - A grim graph object - `grim.Imputation.graph_networkx.Graph`
 - A GRMA donors' graph object - `GRMA.Match.Graph`


```python
from grma.match import matching, find_matches, Graph

donors_graph = Graph.from_pickle("path/to/save.pkl")
matching_results = find_matches("<path to the imputed file>", donors_graph)
```

For more options and information about the input of `matching()` see the documentation in the code.
