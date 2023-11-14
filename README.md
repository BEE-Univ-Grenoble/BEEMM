# BEEMM

## To do list 

### build_markov

*Julie, Matilda, Léo, Laure*

```r
build_markov(sequence_table, order)`
```
- **return** : une table kmer_proba

### markow_likelihood

*Miléna, Lise, Dorine, Manon*

```r 
markow_likelihood(sequence, table_kmer_proba, log=TRUE)`
```

### markov_likelihood_seq

*Lisa, Océane, Justine*

```r 
markov_likelihood_seq(sequence_table, table_kmer_proba)`
```

ou 

```r 
markov_likelihood_seq(sequence_table, ...)`
```

### markov_bayes

```r 
markov_bayes(M0, M1, sequence_testing, priorM0 = 0.5)
```

- Applies markov_likelihood_seq on every M0 and M1
- Computes p(M|Seq) using Bayes
- Returns a tibble similar to  markov_likelihood_seq with four extra columns
  - lprob_M0 and lprob_M1
  - lmodel_M0 and lmodel_M1 

*Julie, Laure, Léo*

### markov_specificity

```r 
markov_specificity(M0s, M1s, sequence_testing, priorM0 = 0.5)
```
- Calls `markov_bayes`
- Computes specificity for all model pairs
- returns 

* Miléna, Lise*
### markov_sensibility

#### function prototype:
```r 
markov_sensibility(M0s, M1s, sequence_testing,  priorM0 = 0.5)
```

#### Computation algorithm
- Computes sensibility for all model pairs by calling `markov_specificity`
  ```
  markov_specificity(M0s = M1s, M1s = M0s, sequence_testing, priorM0 = 1 - priorM0)
  ```

### markov_validate

```r 
markov_validate(sequence_learning_M0, sequence_learning_M1,
                sequence_testing_M0, sequence_testing_M1,
                order_min, order_max,
                priorM0 = 0.5)
```

- prepare the models lists
- call `markov_specificity` & `markov_sensibility`
- merges both results in a three columns table

*Manon, Dorine*

### markov_validate_kfold

```r 
markov_validate_kfold(sequence_learning_M0, sequence_learning_M1,
                      order_min, order_max,
                      priorM0 = 0.5,
                      learning_fraction = 0.8, nrand = 10)
```
- Repeat `nrand` times :
  - splits learning data set in learning and testing data sets
  - calls `markov_validate`
- returns the concatenated results of the `nrand` estimates.

*Lisa, Justine, Oceane*

### viterbi

*Eric*

