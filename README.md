# Data_Science_Course
#Assessment code and data for the data science and machine learning course in Bristol 

# Identification Of Hypothetical Effector Proteins

#Fungal pathogens of plants secrete proteins called effectors to help eade the plant immune system. Effector proteins are generally small and cystine rich. This is because proteins must be small in oredr to fit through thorough small protein channels for secretion out of the cell. Furthermore, as these proteins are secreted outside of the fungal cell, they are suceptible to degredation by proteases. Therefore, a high cysteine content allows the formation of a greater number of disulphide bridges (bonds between cysteines) within the protein's tertiary structure which protects the protein from protease degredation.

#In this report I have generated some code that is able to screen a FASTA file to find proteins with amino acid sequences shorter than 200 amino acids with a cysteine content of between 2 and 20%. This is based on the paper by Lu and Edwards (2016) stating that effector proteins generally have these characteristics.

#Through investigation of effector protein function, we may be able to identify a specific target for control mechanisms for pathogenic fungi which have significant affects in agriculture and healthcare. Finding hypothetical effector proteins through in silico means allows for high throughput screening of entire genomes in a cheap and fast manner. By narrowing down a full species' secretome to a shorter list of potential effectors, we are providing a basis to inform in vivo experiments for effector prediction. This saves, time, money and resources.

# Methods and Results

#In this report I outline some code that can be used to iterate thorugh FASTA files identifying proteins with effector-like properties. Specifically, it creates a new dataframe containing a subset of secreted proteins with sequences shorter than 200 amino acids in length as well as containing a cysteine content of greater than 2%.


#First the libraries needed for the code to run are imported

from Bio import SeqIO

import pandas as pd

import seaborn as sns

import numpy as np

import scipy.stats as stats

from sklearn.model_selection import train_test_split

from sklearn.linear_model import LinearRegression

#and the FASTA file read in 

for record in SeqIO.parse("FastaFile.fa", "fasta"): #fasta file for this is provided through blackboard as data cannot be shared on GitHub
    print(record)

records = list(SeqIO.parse("FastaFile.fa", "fasta"))

print(f"Number of sequences in fasta_file = {len(records)}")

#Here I have created a function that can caluculate the percentage of the amino acid sequence made up of C (cysteine)

def calculate_cysteine_percentage(sequence):
    
    c_count = sequence.upper().count('C')
    
    total_aa_count = len(sequence)
    
    return (c_count / total_aa_count) * 100.0 

#Now I will run through the fasta file and calculate the C content percentage for each sequence as well as its length and put it into a new dataframe called effector_df

#Creating lists to store the cysteine content and length of each sequence

sequence_IDs = []

cysteine_content = []

sequence_lengths = []

#Go through the fasta file using the calculate_cysteine_percentage function

for record in SeqIO.parse("FastaFile.fa", "fasta"):
   
    sequence_ID = record.id
   
    C_content = calculate_cysteine_percentage(record.seq)
    
    sequence_length = len(record.seq) #measures the number of amino acids in sequence

#Adding data into the lists
    
    sequence_IDs.append(sequence_ID)
   
    cysteine_content.append(C_content)
   
    sequence_lengths.append(sequence_length)

#Create a new Data Frame from the lists

effector_df = pd.DataFrame({'Sequence ID': sequence_IDs, 'Cysteine Content': cysteine_content, 'Sequence Length': sequence_lengths})

print(effector_df)

#I will now filter the new effector_df for proteins with a length shorter than 200 amino acids

small_protein = effector_df[effector_df["Sequence Length"] < 200]

#And now combine the length filter with a filter for high cysteine content (over 2%) - these are my predicted effector proteins

predicted_effector = small_protein[small_protein["Cysteine Content"] >2]

predicted_effector

predicted_effector.count()

#This has returned 49 potential effector proteins out of the fasta file



# Data Visualisation

#From the new dataframe of predicted effectors I then decided to look at the distribution of cysteien content and sequence length across the samples.

sns.displot(
   
    data=predicted_effector,
   
    x="Cysteine Content",
)

#Here you can see that most of the predicted effectors have between 2 and 3% cysteine content and there is one sequnce with significantly higher cysteine percentage than the others at 7%. It is possible this is an outlier, however it is still within the 2-20% range specified by Lu and Edwards (2016).

#Next I looked at the distribution of sequence lengths:

sns.displot(
   
    data=predicted_effector,
   
    x="Sequence Length",
)

#This shows a relatively normal distribution of sequence length in the predicted effectors. No sequences were shorter than 60 amino acids. This makes sense as there needs to be a long enough sequence to encode a functioning protein but small enough to allow for secretion.

#My next hypothesis was that there could be a correlation between length and cysteine content in my predicted effectors. It is possible that smaller proteins have higher a cysteine content to avoid digestion in the apoplast. Smaller proteins are more easily digested ad therefore a greater amount of disulphde bridges could help protect small proteins from degreation. To investigate this I plotted a scatterplot to visualise the relationship:

sns.relplot(
   
    data=predicted_effector,
   
    x="Sequence Length",
   
    y="Cysteine Content",
)

#This graph may show a weak negative correlation however it is not very clear.

#I then looked at whether I could create a model to predict the cysteine content based on the sequences length:

X = predicted_effector[["Sequence Length"]] 

y = predicted_effector["Cysteine Content"]

train_X, test_X, train_y, test_y = train_test_split(X, y) #split into test/training set

model = LinearRegression(fit_intercept=True)

model.fit(train_X, train_y)

model.score(test_X, test_y)

#This is a really low model score! It is possible the data wasn't varied enough in length or cysteine content or there weren't enough data points. Alternatively, there may have been no correlation in the first place?

#To see if this was the case I checked the strength of the correlation by calculating the correlation coeficient:

corr = np.corrcoef(predicted_effector["Sequence Length"], predicted_effector["Cysteine Content"])

print(corr)

#This shows a very weak negative correlation (-0.09) between length and cysteine content but this doesnt provide a p-value to show whether this relationship is significant. Instead, below I used the scipy package to calculate the correlation coefficient and the related p-value:

r = stats.pearsonr(predicted_effector["Sequence Length"], predicted_effector["Cysteine Content"])

print(r)

#Now we can see a p-value which shows that this correlation is non-significant (p>0.05) which therefore could explain why the model didn't ahve a very high accuracy score.

# Conclusion 

#The filtering technique used to identify proteins with the main characteristics of aneffector protein and display the %cysteine content as wlel as squence length. From visualisation of the resulting adtaframe, you can see a normal distribution of sequnece lengths under 200 amino acids but the cysteine content was more congregated around the lower percentage content. Using this new dataframe, I created a model to see if cysteine content could be predicted from sequence length based on the hypothesis that smaller proteins would require a higher cysteine content to avoid digestion. However, this model gave a very low score. It is possible that sequence length is in fact not a good predictor of cysteine content however it could also suggest that the dataset was not large enough or did not contain a suitable distribution of sequence lengths and cysteine contents to generate an effective model. While it is possible the latter is true, the extremely low score outputted from the model suggests that it is instead not a strong relationship. For this reaosn, testing the correlation coefficient and resulting p-value provided evidence to suggest that there is no significant correlation between secreted protein length and cysteine content.

#There are a multitude of other features that are common to effector proteins that were not included in this analysis. These could also contribute the the stability of proteins in the apoplast, protecting from digestion by proteases. Therefore, while cystine content is known to be an important factor, other influences could contribute to the correlation being non-significant. Further reserach into ceratin motifs (such as the AXLA motif (Jiang et al., 2016)) that are common to effector proteins could allow for more specific filtering and therefore a more accurate prediction model to be obtained.

#While this code does not provide a perfect prediction method it does significantly narrow down a dataset to a more manageable size for further in vitro investigation. Here 450 amino acid sequences were narrowed down to 49, a significantly more maneagable number.

# References

#Lu, S. and Edwards, M.C., 2016. Genome-wide analysis of small secreted cysteine-rich proteins identifies candidate effector proteins potentially involved in Fusarium graminearum− wheat interactions. Phytopathology, 106(2), pp.166-176.

#Sonah, H., Deshmukh, R.K. and Bélanger, R.R., 2016. Computational prediction of effector proteins in fungi: opportunities and challenges. Frontiers in plant science, 7, p.126.
