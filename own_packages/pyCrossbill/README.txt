###################################################################################################################

Hi I'm Aniket yadav,

	Welcome to my python package which is a type of computational biology and bioinformatics tool, this is a program based
	package not includes user interface now but, I can try with UI also at next time.
		According to NCBI, Bioinformatics is defined as the application of tools of computation and analysis to 
	the capture and interpretation of biological data. It is an interdisciplinary field, which harnesses computer 
	science, mathematics, physics, and biology. Bioinformatics is essential for management of data in modern biology 
	and medicine.
	This is the easy and suitable way to analyse a large number of DNAs, RNAs and protein sequence.
	Package contains number of folder like:
		
		|-- pyCrossbill
			|-- BioLearning
				|-- __init__.py
				|-- Learning.py
			|-- Biological_Tool
				|-- __init__.py
				|-- DistanceTools.py
				|-- GeneralChem.py
				|-- qualitySequence.py
				|-- SequenceTools.py
			|-- database
				|-- __init__.py
				|-- Compound_reterive_.py
				|-- DBTool.py
				|-- pyCrossbill.sqlite  (default Sqlite Database File)
			|-- datasets
				|-- drugset.csv
				|-- raw_drug.csv
			|-- docs
				|-- samples
					|-- cov2-Wuhan-Hu1.fasta
					|-- dengu3_strain-D00-0107.fasta
					|-- ebola-HomoSapiens_India-.fasta
					|-- nipah-G-gglycoproteinImmunogenic.fasta
			|-- Formats
				|-- __init__.py
				|-- FileFormat.py
			|-- __init__.py
			|-- sys.py
			|-- License.txt
			|-- README.txt
			|-- setup.py

By default package is use sqlite database to store local data of sequence or use to make your own biological database
this is the package contains so many types of methods and classes to do a simple task, 

-----print(load_csv(file).dataset())

output: [['Aspirin', 'CC(=O)OC1=CC=CC=C1C(=O)O', '135'], ['Sevoflurane', 'C(OC(C(F)(F)F)C(F)(F)F)F', '25'], 
['Propofolï¿½', 'CC(C)C1=C(C(=CC=C1)C(C)C)O', '18'], ['Diazepam', 'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3', '132'], 
['Midazolam', 'CC1=NC=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4F', '159'], 
['Morphine Sulfate', 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O.OS(=O)(=O)O', '255'], 
['Promethazine', 'CC(CN1C2=CC=CC=C2SC3=CC=CC=C31)N(C)C', '60'], ['Carbamazepine', 'C1=CC=C2C(=C1)C=CC3=CC=CC=C3N2C(=O)N', '190.2'], 
['Diazepam', 'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3', '132'], ['Magnesium sulfate', '[O-]S(=O)(=O)[O-].[Mg+2]', '1124'], 
['Phenobarbitone', 'CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2', '174'], 
['Phenytoin Sodium', 'C1=CC=C(C=C1)C2(C(=O)[N-]C(=O)N2)C3=CC=CC=C3.[Na+]', 'n a.'], 
['Sodium valproate', 'CCCC(CCC)C(=O)[O-].[Na+]', 'n.a.'], 
['lorazepam', 'C1=CC=C(C(=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O)Cl', '167']]

-----print(load_csv(file).dataset(column=1, upto=5))

output: ['CC(=O)OC1=CC=CC=C1C(=O)O', 'C(OC(C(F)(F)F)C(F)(F)F)F', 'CC(C)C1=C(C(=CC=C1)C(C)C)O', 
'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3', 'CC1=NC=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4F']

class statistics contains many types of methods which is the statistics functions like, (mean, median, variance, standard deviation)
we can implements those programs as:
 
-----set = [3, 6, 9, 2, 7]
-----print(statistics(set).standard_dev())

output: 2.5768197453450252

-----print(statistics(set).mean())

output: 5.4

we need a way to determine if thhere is linear correlation or not, so we calculate what is know  
as the "PRODUCT-MOMENT CORRELATION COEFFICIENT".

-----pd = [12.3, 53.2, 52.2, 13.4, 83.5]
-----at = [26.7, 34.2, 52.0, 63.5, 12.6]
-----print(Product_moment_CC(at, pd))

output: -402.4764650551727

-----pd = [12.3, 53.2, 52.2, 13.4, 83.5]
-----at = [26.7, 34.2, 52.0, 63.5, 12.6]
-----print(covarince(at, pd))

output: -269.32599999999996


After this, we'll move toward distances between two sequences, we define this distance as the number
of bases by which they differ.

firstly, find lemple Ziv between two sequence:
return two outputs unique list of codons/elements and and total number of element in list.

-----seq1 = 'CATGTG'
-----seq2 = 'CATGTT'
-----both = 'CATGTGCATGTT'
-----seq2 = LempeleZiv(seq2)
-----both = LempeleZiv(both)
-----seq1 = LempeleZiv(seq1)

output: (['C', 'A', 'T', 'G', 'TT'], 5)

lemple Ziv is commonly use in find Normalize compression distance so now we finded distance like this:

-----print(NormalizeCompressionDistance(seq1, seq2, both))

output: 0.4
*both argument is defines the merge up of both seq1 and seq2 sequences

finding euclidean distance between two sequence

-----print(EuclideanDistance(seq1, seq2))

output: 1.4142135623730951

