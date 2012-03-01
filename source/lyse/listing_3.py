from lyse import *
from pylab import *

# Let's obtain the dataframe for all of lyse's currently loaded shots:
df = data()

# Now let's see how the MOT load rate varies with, say a global called
# 'detuning', which might be the detuning of the MOT beams:

detunings = df['detuning']

# mot load rate was saved by a routine called calculate_load_rate:

load_rates = df['calculate_load_rate', 'mot loadrate']

# Let's plot them against each other:

plot(detunings, load_rates,'bo',label='data')

# Maybe we expect a linear relationship over the range we've got:
m, c = linear_fit(detunings, load_rates)
# (note, not a function provided by lyse, though I'm sure we'll have
# lots of stock functions like this available for import!)

plot(detunings, m*detunings + c, 'ro', label='linear fit')
legend()

#To save this result to the output hdf5 file, we have to instantiate a
#Sequence object:
seq = Sequence(path, df)
seq.save_result('detuning_loadrate_slope',c)
