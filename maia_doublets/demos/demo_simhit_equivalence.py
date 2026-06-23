import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

f_lcio = "background_sim_lcio/simhits.pkl"
f_root = "background_sim_root/simhits.pkl"

lcio_hits = pd.read_pickle(f_lcio)
root_hits = pd.read_pickle(f_root)

print(root_hits.columns)
print(lcio_hits.columns)
print("np.all(root_hits.columns == lcio_hits.columns)\n",
      np.all(root_hits.columns == lcio_hits.columns))

try:
    assert_frame_equal(root_hits, lcio_hits, check_dtype=False, check_column_type=False)
    print("Hits dataframes are equal according to assert_frame_equal")
except AssertionError as e:
    print(f"Hits dataframes are not equal: {e}")

