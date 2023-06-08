# created a class for reading FASTQ file format
# with takes an argument of file location from the local machine
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

    # function to return length of the FASTQ file 
    # read sequences in single file
    def FASTQlen(self):
        fastq = self.readFASTQfile()   
        # returns the length of the dict...
        return len(fastq)
    
    # function for returning top most limit of sequence where default value of top reads is 5
    def head(self, top=5):
        # takes an argument 'top' for reading or returning limited sequence from top 
        fastq_limited = self.readFASTQfile()
        # read and select all fastq and makes a list of total tuple items
        lis_limit = list(fastq_limited.items()) 
        return dict(lis_limit[:top])  # and create a dictionary of top limited list
 

f = 'A:/PROJECTS/out_1000_ERR101899.1.fastq'
fastq_ = 'A:/PROJECTS/ERR101900_1.fastq'
FASTQ = FASTQformat(fastq_)

# print(FASTQ.readFASTQfile())
# FASTQlen() returns 0 if file will empty..
print(FASTQ.FASTQlen())
# head(top=3) -> return three reads from top [return {} if file will empty]
print(FASTQ.head(top=3))