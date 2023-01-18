# grma
> grma - GRaph based MAtching

grma is a package for finding HLA matches using graphs approach.
The matching is based on [grim's](https://github.com/nmdp-bioinformatics/py-graph-imputation) imputation.

You should follow these steps for finding matches:
* Build 'Donors Graph' - your HLA search space object.
* Impute your patient's (or any other) genotypes with grim algorithm
(You may use the grim's default settings or use your own settings).
* Use grma algorthm for finding matches efficiently.

### Building The Donors' Graph:

The donors' graph is a graph which contains all the donors (the search space). It implemented using a LOL (List of Lists) representation written in cython for better
time and memory efficiency.
The building might take a lot of memory and time, so it's recommended to save the graph in a pickle file.

Before building the donors' graph, all the donors' HLAs must be imputed using `grim`.
Then all the imputation files must be saved under the same directory.

```python
from grma.donorsgraph.build_donors_graph import BuildMatchingGraph

PATH_TO_DONORS_DIR = "./data/donors_dir"
PATH_TO_DONORS_GRAPH = "./data/donors_graph.pkl"

build_matching = BuildMatchingGraph(PATH_TO_DONORS_DIR)
graph = build_matching.graph  # access the donors' graph

build_matching.to_pickle(PATH_TO_DONORS_GRAPH)  # save the donors' graph to pickle
```

### Imputing patients' genotypes:
The function `matching` apply both grim and grma algorithms.
It gets a path to a grim configuration file with the settings of the algorithm and the path to the data files.
You can't ignore this configuration file, then grim will work with default settings (Notice: If you will use the default configs, you will have to put your data files in a specific path according to grim's default config).

If you already have imputed genotypes, you can use the function `find_matches` which will apply only the grma algorithm.

### Search & Match

The functions `matching` \ `find_mathces` find matches up to 3 mismatches and return a `pandas.DataFrame` object of the matches sorted by number of mismatches and their score.

They get these parameters:
* imputation_filename: a path to the file of the patients' typing. (only in `find_matches`)
* grim_config_file: a path to `grim` configuration file (optional, only in `matching`)
* match_graph: a grma donors' graph object - `grma.match.Graph`

```python
from grma.match import Graph, find_matches

PATH_TO_PATIENTS_FILE = "./data/patients_file.txt"
PATH_TO_DONORS_GRAPH = "./data/donors_graph.pkl"

# The donors' graph we built earlier
donors_graph = Graph.from_pickle(PATH_TO_DONORS_GRAPH)
matching_results = find_matches(PATH_TO_PATIENTS_FILE, donors_graph)

# matching_results is a dict - {patient_id: the patient's result dataframe}

for patient, df in matching_results.items():
    # Use here the dataframe 'df' with the results for 'patient'
    print(patient, df)
```

`find_matches` and `matching` take some optional parameters, which you might want to change:

* search_id: An integer identification of the search. default is 0.
* donors_info: An iterable of fields from the database to include in the results. default is None.
* threshold: Minimal score value for a valid match. default is 0.1.
* cutof: Maximum number of matches to return. default is 50.
* verbose: A boolean flag for whether to print the documentation. default is False
* save_to_csv: A boolean flag for whether to save the matching results into a csv file. default is False.
* calculate_time: A boolean flag for whether to return the matching time for patient. default is False.
  In case `calculate_time=True` the output will be dict like this: `{patient_id: (results_dataframe, time)}`


### Set Database
In order to get in the matching results more information about the donors than the matching information,
one can set a database that has all the donors' information in it.
The database must be a `pandas.DataFrame` that its indexes are the donors' IDs.

After setting the database, when calling one of the matching functions,
you may set in the `donor_info` variable a `list` with the names of the columns you want to join to the result dataframe from the database.

Example of setting the database:
```python
import pandas as pd
from grma.match import set_Database

donors = [0, 1, 2]
database = pd.DataFrame([[30], [32], [25]], columns=["Age"], index=donors)

set_Database(database)
```

Note: `set_Database()` need to be called only once before
performing the matching if additional fields are ment to be added to the results.
