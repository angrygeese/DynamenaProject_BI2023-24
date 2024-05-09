import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from IPython.display import display


pd.options.display.float_format = '{:,.2f}'.format


def run_eda(df: pd.DataFrame, factor_thershold: int=35, do_na_count=True):

    print('Hi there, fellow researcher.', end='\n\n')

    df_shape = df.shape
    print(f'Dataframe consists of {df_shape[1]} variables and {df_shape[0]} observations', end='\n\n')

    col_data_types = df.dtypes.to_dict()
    data_types = {'string': [], 'number': [], 'factor': []}
    for key, value in col_data_types.items():
        if df[key].nunique() <= factor_thershold:
            data_types['factor'].append(key)
        else:
            if value == np.object_:
                data_types['string'].append(key)
            elif value == np.int64 or value == np.float64:
                data_types['number'].append(key)
    for var_type, var in data_types.items():
        print(f'Following variables can be classified as {var_type}-type: {*var,}')
    print(f'Criteria are unique values count (lower than {factor_thershold} is factor) and then data type ("int64"/"float64" is number, "object" is string).', end='\n\n')

    fig, axes = plt.subplots(nrows=2, ncols=len(data_types['number']), figsize=(24, 6), dpi=300)

    print('1) Basic statistics for factor-type variables.', end='\n\n')
    if data_types['factor']:
        for column in data_types['factor']:
            col = df[column]
            col_value_counts = col.value_counts()
            col_frequencies = col_value_counts / col.size * 100
            print(f'For variable "{column}:\n- Counts:\n{col_value_counts.to_string()}\n- Frequencies:\n{col_frequencies.to_string()}%', end='\n\n')
    else:
        print('No factor-type variables.', end='\n\n')

    print('2) Basic statistics for numeric variables.')
    
    num_stats = df[data_types['number']].describe()
    col_outliers, col_n_unique = [], []

    for pos, column in enumerate(data_types['number']):
        col = df[column]
        col_mean, col_IQR = num_stats.at['mean', column], num_stats.at['75%', column] - num_stats.at['25%', column]
        col_outliers.append(col[(col < col_mean - 1.5 * col_IQR) | (col > col_mean + 1.5 * col_IQR)].size)
        col_n_unique.append(col.nunique())
        sns.boxplot(x=col, ax=axes[0, pos])
        sns.histplot(col, bins='auto', ax=axes[1, pos])
        
    num_stats.loc['n_outliers'] = col_outliers
    num_stats.loc['n_unique'] = col_n_unique

    display(num_stats)

    print('Numeric variables visualization:', end='\n')
    
    plt.show()

    if do_na_count:
        
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(18, 4), dpi=300)

        df_na = df.isna()
        na_count = df_na.sum().sum()
        na_var_count = dict(zip(df.columns, df_na.sum(axis=0)))
        na_rows = df_na.any(axis=1).sum()
        na_columns = df.columns[df_na.any(axis=0)]
        print('3) NA values:')
        print(f'- Total NA count is {na_count} in {na_rows} rows.\n- Following columns contain NA: {*na_columns, }', end='\n\n')
        sns.barplot(pd.DataFrame(na_var_count, index=['NA_count'], columns=pd.Series(df.columns, name='Variables')), ax=axes[0])
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=90)
        axes[0].set_title("NA count")
        None

        duplicates_count = df.duplicated().sum()
        print(f'Dataframe contains {duplicates_count} duplicates.')


        sns.heatmap(df[data_types['number']].corr(), cmap='mako', ax=axes[1])
        axes[1].set_title("Correlation matrix for numerical variables")
        None

        plt.show()