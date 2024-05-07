# GeneClassCoTrain
Term project for Big Data Science (CSCI 5952) - Pilot project: Using a cotraining semi-supervised method for gene classification

## **1. propagate_hpo_genes.R**
**Input**: HPO data
**Output**: Propagated gens associated to phentype of diseases of interest. For clas project it creates a TSV file contianing genes associated to 
TB-related phenotypes. 

**Files needed**: 
1. hp.obo
2. mondo_2024_03.obo
3. phenotype.hpoa
4. genes_to_phenotype.txt

**Functions**
1. map_mondo_to_ontology (HPO df, obo file): Maps MONDO disease identifiers to the HPO disease identifiers (OrphaNet and OMIM).
2. get_ances (list of disease identifiers, obo file): gets the ancestor term of all the terms found in the hpo data set.
3. propagate_genes (list of disease identifiers from HPO data, ancester df, df containing disease id and their annotated hp ids, df with mapped MONDO identifiers, MONDO obo df): function that propagates genes from children terms up to parent terms.

**Output**: 
1. propagated_hpo_disgenes.RData
2. propagated_hpo_genes.tsv

## **2. propagate_disgenet_genes.R**
**Input**: Disease gene association data form DisGeNET and MONO obo file.
**Output**: Propagated genes associations for the disease of interest (TB).

**Files needed**: 
1. mondo_2024_03.obo
2. all_gene_disease_associations.tsv

**Functions**
1. map_mondo_to_ontology (DisGeNET df, obo_file): Maps MONDO disease identifiers to the disease identifiers contained in DisGeNET database(UMLS).
2. get_ances (list of disease identifiers, obo file): gets the ancestor term of all the terms found in the DisGeNET data set.
3. propagate_genes (list of disease identifiers from DisGeNET df, ancester df, DisGeNET data with mapped MONDO identifiers): function that propagates genes from children terms up to parent terms.

**Output**: 
1. tb_propagated_disgenet_genes.tsv
2. all_disgenet_prop_genes.RData

## **3. final_term_proj_model.ipynb**
**Input**: BioGRID data, propagated genes obtained from DisGeNET and HPO.
**Output**: text file containing performance results for the two cotraining models. 

**Files needed**: 
1. propagated_hpo_genes.tsv
2. tb_propagated_disgenet_genes.tsv
3. biogrid_network.txt

**Functions**
1. get_embeddings (adjacency matrix): creates embeddings for the adjacency matrix.
2. genes_in_netowork(positive gene list, adjacency matrix genes to index, adjacency matrix embeddings): function checks to see if the positive genes are in the network. If they are not then thye are taken out. Takes the rows corresponding to the positive list to create embeddings contining positive examples. A dictionary containing where the genes are indexed is outputed with the positive embeddings. 
3. get_negatives(positive gene list, djacency matrix genes to index, adjacency matrix embeddings): like the genes_in_network but instead embeddings for negatives is obtained. Only genes in network are used as negatives genes and negative genes are picked randomly. Same number of negatives are obtained as there are positives. A dictionary containing where the genes are indexed is outputed with the negative embeddings. 
4. get_feat_matrix (embeddings containing positive genes, embeddings containing negative genes, gene to index for the positive embeddings, gene to index for the negative embeddings): joins th enegative and positive embeddings to create the feature matrix for training. A dictionary containing where the genes are indexed is outputed with the feature matrix embeddings. 
5. get_target_vector (feature embeddings, gent to index dictionary of the feature embeddings): creates the target vector for the feature matrix.
6. logistic_regression(): parameter settings for the logistic regression classifier. 
7. random_forest(): parameter settings for the random forest classifier.
8. label_propagation(): parameter settings for the label propagation classifier.
9. co_train_lr_rf (logistic regression object, random forest object, disgenet training split, hpo training split, disgenet training labels, hpo labels, disgenet unlabeld data, hpo unlabeled data, disgenet test split, hpo test split, disgenet test split labels, hpo test labels, threshold confidence setting, adjacency matrix embeddings): creates cotrain 1 model instance. Training and testing phase using DisGeNET and HPO unlabeled and labeled data. Performance metreics are calculated and placed in dictionary. 
10. co_train_lr_label_prop(logistic regression object, label propagation object, disgenet training split, hpo training split, disgenet training labels, hpo labels, disgenet unlabeld data, hpo unlabeled data, disgenet test split, hpo test split, disgenet test split labels, hpo test labels, threshold confidence setting, adjacency matrix embeddings): creates cotrain 2 model instance. Training and testing phase using DisGeNET and HPO unlabeled and labeled data. Performance metreics are calculated and placed in dictionary. 
