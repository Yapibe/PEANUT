import pandas as pd
from os import path


def shuffle_scores(dataframe, shuffle_column):
    dataframe[shuffle_column] = dataframe[shuffle_column].sample(frac=1).reset_index(drop=True)
    return dataframe


def main():
    # Load the scores from the Excel file
    main_dir = path.dirname(path.dirname(path.abspath(__file__)))
    scores_file = path.join(main_dir, 'Inputs', 'experiments_data', 'scores_T_v_N.xlsx')
    scores_df = pd.read_excel(scores_file)

    # Shuffle the scores
    # Assuming 'Score' is the column name; replace with the actual column name if different
    shuffled_scores_df = shuffle_scores(scores_df, 'Score')

    # Save the shuffled scores to a new Excel file
    shuffled_scores_df.to_excel('shuffled_scores.xlsx', index=False)


if __name__ == "__main__":
    main()
