import numpy as np
import pandas as pd



df = pd.read_csv("Files/dndeta.csv", names=['Participates', 'Values'])

weights = np.repeat([1.0, 5.0], [5, 15])

print(f"Weights: {len(weights)}")
print(f"Values: {len(df['Values'])}")
print(f"Participates: {len(df['Participates'])}")

sum_value = np.sum(df["Participates"] * df["Values"] * weights)
sum_weights = np.sum(weights)
sum_participates = np.sum(df["Participates"] * weights)

average = sum_value / sum_weights
average_part = sum_participates / 80

print(sum_participates)
print(f"Average: {average/12.47} with {average_part} Participates")

