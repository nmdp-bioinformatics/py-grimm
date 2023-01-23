from grma.donorsgraph.build_donors_graph import BuildMatchingGraph
from grma.match import find_matches

DONORS_DIR_PATH = "./data/test_donors"
PATIENTS_FILE_PATH = "data/test_patients.txt"

build_matching = BuildMatchingGraph(DONORS_DIR_PATH)
donors_graph = build_matching.graph

matching_results = find_matches(
    PATIENTS_FILE_PATH,
    donors_graph,
    threshold=0.2,
    cutof=100,
    calculate_time=False,
    verbose=False,
)

for patient, df in matching_results.items():
    # There are only few donors and patients (The donors\patients files are the same),
    # so the best match for the patient should be itself.
    assert df["Donor_ID"][0] == df["Patient_ID"][0]
