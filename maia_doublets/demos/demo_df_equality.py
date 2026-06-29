"""
This script checks if two pandas DataFrames are equal by loading them from pickle files.
It provides an option to sort the DataFrames by a specified column before comparison.
If the DataFrames are equal, it prints a confirmation message; otherwise, it prints the differences.
"""
import argparse
import pandas as pd
from pandas.testing import assert_frame_equal
from pandas.testing import assert_series_equal

def main():
    args = options()
    df1 = pd.read_pickle(args.df1)
    df2 = pd.read_pickle(args.df2)
    print(f"Size of df1: {df1.shape}, Size of df2: {df2.shape}")
    if args.sortby:
        sortby = args.sortby.split(",")
        df1 = df1.sort_values(by=sortby).reset_index(drop=True)
        df2 = df2.sort_values(by=sortby).reset_index(drop=True)
    try:
        assert_frame_equal(df1, df2, check_dtype=False, check_column_type=False, check_like=True)
        print("DataFrames are equal.")
    except AssertionError as e:
        print("DataFrames are not equal.")
        print(e)
        # check for differences column-by-column
        for col in sorted(df1.columns):
            if col in df2.columns:
                try:
                    assert_series_equal(df1[col], df2[col], check_dtype=False, check_names=False)
                except AssertionError as e:
                    print(f"Column '{col}' is different:")
                    # print(e)
            else:
                print(f"Column '{col}' is missing in the second DataFrame.")


def options():
    parser = argparse.ArgumentParser(description="Check if two dataframes are equal.")
    parser.add_argument("df1", type=str, help="Path to the first dataframe pickle file.")
    parser.add_argument("df2", type=str, help="Path to the second dataframe pickle file.")
    parser.add_argument("--sortby", type=str, default=None, help="Optional column name to sort the dataframes before comparison.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
