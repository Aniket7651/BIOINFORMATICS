
# created a class for reading FASTQ file format
# with takes a argument of file location on the local machine
class FASTQformat:

    def __init__(self, file_Location):
        # initialize argument of file location
        self.fileLocation = file_Location

    # function of reading FASTQ from FASTQformat class, returns dictionary of sequence and thier respective ASCII characters
    # ie. dict = {'seq': 'ASCII', ....}
    def readFASTQfile(self):
        # seqData is the empty list use for storing data of sequence and ASCII only,
        # seqData variable drop all the sep (+) and header line starting from '@'
        seqData = []
        # seq and ascii list variable stores sequence and thier ASCII code seperately...
        (seq, ascii) = ([], [])
        with open(self.fileLocation, 'r') as FASTQ:
            # open FASTQ file and read as 'FASTQ', reading start from 2nd lines and store in 'line' variable
            # ie. try to drop first header.
            line = FASTQ.readlines()[1:]
            for i in range(len(line)):
                if i%2 == 0:         # seqData stores those line which contain ASCII and sequence both (even lines)
                    seqData.append(line[i].replace('\n', ''))    # replace all new lines (\n) of the each lines
        
        # this is the seperate part of this program, makes a 'seq' and 'ascii' list from 'seqData', and than 
        # a dict (seqdict) from 'seq' and 'ascii'...
        for i in range(len(seqData)):
            if i%2 == 0:       # this part is same as above even line part, where 'seq' store sequence data from 'seqData'
                seq.append(seqData[i])
            elif i%2 != 0:      # except all (ASCII CODE) stores 'ascii'
                ascii.append(seqData[i])
        seqdict = {key: value for key, value in zip(seq, ascii)}  # use dict. comprehension to make dict of 'seq' and 'ascii'
        # return dict of final sequence and thier respective ASCII code for finding the phred score of that sequence
        return seqdict

    # fuction to return length of the FASTQ file 
    # read sequences in single file
    def FASTQlen(self):
        fastq = self.readFASTQfile()   
        # reutrns the length of the dict...
        return len(fastq)
    
    # function for returning top most limit of sequence
    def head(self, top=5):
        # takes a argument 'top' for reading or returning limited sequence from top 
        fastq_limited = self.readFASTQfile()
        # read and select all fastq and makes a list of total tuple items
        lis_limit = list(fastq_limited.items()) 
        return dict(lis_limit[:top])  # and create a dictionary of top limited list


def FASTQ_to_FASTA_Convertion(inputFile, outputFile, limits=5):
    import re
    with open(inputFile, 'r') as readFastsq:
        readlist = readFastsq.readlines()
        (header, seqs) = ([], [])
        for line in range(len(readlist)):
            if line%2 == 0:
                if readlist[line][0] == '@':
                    replacesign = readlist[line].replace('@', '>')
                    header.append(replacesign)
            else:
                if re.search("^[ATGC]+$", readlist[line]): 
                    seqs.append(readlist[line])
    
    outfile = f"{outputFile}out_{header[0].split()[0][1:]}.fasta"
    with open(outfile, 'a') as outfile:
        l = 0
        for h, s in zip(header, seqs):
            if l < limits:
                outfile.write(h)
                outfile.write(s)
                l += 1
            else: break
        outfile.close()
    return outfile

# print(FASTQ_to_FASTA_Convertion('A:/PROJECTS/BIOINFORMATICS_software/own_packages/pyCrossbill/test.txt', 'A:/PROJECTS/BIOINFORMATICS_software/own_packages/pyCrossbill/'))

# there is the seperate function for removing duplicates from the list of the sequences
# returns the tuple of index and unique
def remove_duplicates(sequence_list):
    # 'index', where only unique sequence is counted and 'unique', filter the item in list of sequence_list variable
    (unique, index) = ([], [])
    for i in sequence_list:
        if sequence_list.index(i) not in index:
            index.append(sequence_list.index(i))
    
    for i in index:
        unique.append(sequence_list[i])

    return index, unique


def GC_percentage(seq_dict):
    GC_ = []
    def calc_GC_percentage(seq):
        (count_g, count_c) = (0, 0)
        for code in seq:
            if code == 'G':
                count_g += 1
            if code == 'C':
                count_c += 1
        return 100 * float((count_g + count_c)) / len(seq)

    for key in seq_dict:
        GC_.append(round(calc_GC_percentage(key), 2))
    return GC_


def transform(mtrix):
    transf = []
    for i in range(len(mtrix[0])):
        lisr = []
        for j in range(len(mtrix)):
            lisr.append(mtrix[j][i])
        transf.append(lisr)
    return transf


class Quality_scoring:

    def __init__(self, fastq_seq):
        self.seq_ = fastq_seq
    
    def pred_score_33(self, str_values):
        pred_list = []
        for char in str_values:
            pred_list.append(ord(char)-33)
        return pred_list

    def pred_score_64(self, str_values):
        pred_list = []
        for char in str_values:
            pred_list.append(ord(char)-64)
        return pred_list

    def isBASE(self):
        base_ = 'ILLUMINA_33'
        for k in self.seq_:
            ascii = self.seq_[k]
            for char in ascii:
                if char.islower() == True:
                    base_ = 'ILLUMINA_64'
                break
        return base_

    def quality(self):
        reading_scores = []
        for k in self.seq_:
            
            if self.isBASE() == 'ILLUMINA_33':
                reading_scores.append(self.pred_score_33(self.seq_[k]))
            if self.isBASE() == 'ILLUMINA_64':
                reading_scores.append(self.pred_score_64(self.seq_[k]))
        trasf = transform(reading_scores)
        mean = []
        for item in trasf:
            mean.append(sum(item)/len(item))
        return mean, trasf

    def best_readedScore(self):
        quality_list = Quality_scoring(self.seq_).quality()
        majorError = []
        for phreds in quality_list:
            error_prob = []
            for phred in phreds:
                P = 10**(-phred/10)
                error_prob.append(round(P, 6))
            S = sum(error_prob)/len(error_prob)
            majorError.append(round(S, 6))
        return majorError.index(min(majorError)), quality_list[majorError.index(min(majorError))]


class visualization:

    def GC_graph(self, gc_lis, type='scatter'):
        import matplotlib.pyplot as plt
        import seaborn as sns
    
        Y = [i for i in range(len(gc_lis))]
       
        if type == 'dist' or type == 'hist':
            plt.title('GC CONTENT PER SEQUENCE DISTRIBUTION', font='arial', color='#636057', fontsize=18)
            sns.distplot(gc_lis, kde=True, color='#f78181', bins=30)
            
        else:
            plt.title('GC% PER SEQUENCE READS', font='arial', color='#636057', fontsize=18)
            plt.scatter(gc_lis, Y, color='black', facecolor='#f78181')
            plt.xlabel('GC in %', font='arial', color='#636057', fontsize=12)
            plt.xlim(0, 100)
            plt.ylabel('Number of reads', font='arial', color='#636057', fontsize=12)
        return plt.show()
    
    def scoring_graph_BASE33(self, mean, data, style='default'):
        import matplotlib.pyplot as plt
        import seaborn as sns

        nt = [i for i in range(len(mean))]
        if style == 'gray':
            # (span1 0-20,  span2 20-28,  span3 28-42,  mean line color,  box color)
            colors = ('#999999', '#cccccc', '#f2f2f2', '#000000', 'white')
        elif style == 'white':
            colors = ('white', 'white', 'white', 'black', '#bfbdbd')
        elif style == 'cool':
            colors = ('#546af7', '#5ca1fa', '#add6f7', '#112ff0', '#f2f207')
        elif style == 'hot':
            colors = ('#f2b750', '#f7cb7e', '#fce4bb', '#f70505', '#f79d8d')
        elif style == 'heatmap':
            colors = ('#65c9f7', '#9ef6f7', '#f79e9e', '#0f0f0f', '#f2a1f7')
        else: 
            colors = ('#ed7272', '#ecfa82', '#82fa8c', '#ff5c33', 'yellow')

        ax = plt.subplots()[1]
        plt.title(f'SCORING GRAPH (ASCII BASED 33), length {len(mean)} bp', size=20, font='arial', color='#636057')
        plt.plot(nt, mean, c=colors[3], linewidth=1)
        sns.boxplot(data, showfliers=False, width=0.9, color=colors[4], linewidth=0.8)
        ranges = [i for i in range(0, len(mean), 5)]
        ax.grid(1)
        ax.margins(0)
        ax.axhspan(0, 20, facecolor=colors[0], alpha=0.5)
        ax.axhspan(20, 28, facecolor=colors[1], alpha=0.5)
        ax.axhspan(28, 42, facecolor=colors[2], alpha=0.5)
        ax.set_xticks(ranges)
        ax.set_xticklabels(ranges)
        plt.xlim(0, len(nt))
        plt.ylim(0, 42)

        plt.xlabel('Nucleotides (bp)', font='arial', fontsize=12, color='#636057')
        plt.ylabel('Phred Score (Q)', font='arial', fontsize=12, color='#636057')
        return plt.show()


class Trimming:
    
    def __init__(self, FASTQ_file):
        self.fastq = FASTQ_file

    def trimFASTQ_reads(self, chunks, outputfileLocation='null'):
        with open(self.fastq, 'r') as FASTQ:
            line = FASTQ.readlines()
            seqlines = line[:chunks*4]
        if outputfileLocation == 'null':
            file = self.fastq.split('/')
            file.pop()
            file.append(f'out_{chunks}_{seqlines[0].split()[0][1:]}.fastq')
            output = '/'.join(file)
        else: output = f'{outputfileLocation}out_{chunks}_{seqlines[0].split()[0][1:]}.fastq'

        with open(output, 'a') as outputFASTQ:
            for li in seqlines:
                outputFASTQ.write(li)
        return output
    
    def seqTrimming(self, probs, outputLocation):
        lines = ''
        with open(self.fastq, 'r') as FastqR:
            r = FastqR.readlines()
            for line in range(len(r)):
                if line%2 != 0:
                    even = r[line]
                    if '-' in probs:
                        frto = probs.split('-')
                        lines += even[int(frto[0]):int(frto[1])]+'\n'
                    else: lines += even[:probs]+'\n'
                else: lines += r[line][:-1]+f' => length={probs}\n'
                
        outfile = f"{outputLocation}out_{probs}trim{r[0].split()[0]}.fastq"
        with open(outfile, 'a') as fastqOut:
            fastqOut.write(lines)
            fastqOut.close()
        return outfile


########################################## END OF THE PROGRAMS ###############################################

file_path = 'A:/PROJECTS/BIOINFORMATICS_software/own_packages/pyCrossbill/test.txt'
fastq_ = 'A:/PROJECTS/ERR101899.1_sample.fastq'
# FASTQ = FASTQformat('A:/PROJECTS/out_140trim@ERR101899.1.fastq')
FASTQ = FASTQformat('A:/PROJECTS/ERR101900_1.fastq') # 
s1 = 'CCTACGGGTGGCAGCAGTGAGGAATATTGGTCAATGGACGGAAGTCTGAACCAGCCAAGTAGCGTGCAG'
s2 = 'CCTACGGGTGGCTGCAGTGAGGAATATTGGACAATGGTCGGAAGACTGATCCAGCCATGCCGCGTGCAG'
# print(remove_duplicates(['ACGTGGCTGATCGAT', 'ACATGCGGGATCGAT', 'ACGTGGCTGATCGAT', 'AAGGATCATCGATTT', 'ACATGCGGGATCGAT']))
lis = ['ACGTGGCTGATCGAT', 'ACATGCGGGATCGAT', 'ACGTGGCTGATCGAT', 'AAGGATCATCGATTT', 'ACATGCGGGATCGAT']

# dict = FASTQ.head(500)
# quality = Quality_scoring(dict).quality()
# print(quality)
# visualization().scoring_graph_BASE33(quality[0], quality[1], 'cool') # style=['hot', 'cool', 'heatmap', 'white', 'gray']
# print(Trimming(fastq_).trimFASTQ_reads(1000, 'A:/PROJECTS/BIOINFORMATICS_software/own_packages/pyCrossbill/'))
# print('number of the readings: ', FASTQ.FASTQlen())
# gc = GC_percentage(dict)
# print(len(gc))
# print(Trimming('A:/PROJECTS/out_1000_ERR101899.1.fastq').seqTrimming('14-140', 'A:/PROJECTS/'))
# visualization().GC_graph(gc)