"""
Name: Richard Ahn
CS 295 Winter 2012
HW6.py
"""

import re



class ReadMapper():
    """A ReadMapper is a class that handles finding the locations of short
    sequences of DNA (reads) in a very long sequence of DNA (reference).
    
    For simplicity, a seed-and-extend approach is taken, meaning a portion of
    the read must match exactly to the reference genome.
    
    """
    def __init__(self, reference, seed_length, method='hamming'):
        """Create a new ReadMapper, with a fixed reference and seed_length.
        
        You should use build_index to build an index on this reference and call
        the resulting dictionary `self.ref_index`
        
        You should also store the seed_length as `self.seed_length`
        
        Finally, you should store the reference sequence as `self.reference`
        
        """
        self.reference = reference
        self.seed_length = seed_length
        self.ref_index = {}

    def sliding_window(self, reference, window_length):
        """Create a sliding window on reference, returning a list of all substrings
        in the reference.
        
        You should convert the reference to upper-case characters before
        working on it. You should NOT REPORT any substrings that have an N character in them.        
        
        For example:
        >>> reference = 'AAATTTGGGCCC'
        >>> mapper = ReadMapper(reference, 3)
        >>> mapper.sliding_window(reference, 3)
        ['AAA', 'AAT', 'ATT', 'TTT', 'TTG', 'TGG', 'GGG', 'GGC', 'GCC', 'CCC']
        >>> mapper.sliding_window(reference, 5)
        ['AAATT', 'AATTT', 'ATTTG', 'TTTGG', 'TTGGG', 'TGGGC', 'GGGCC', 'GGCCC']
        
        # make sure you don't report the N character
        >>> reference = 'AAATTNGGGCCC'
        >>> mapper = ReadMapper(reference, 3)
        >>> mapper.sliding_window(reference, 3)
        ['AAA', 'AAT', 'ATT', 'GGG', 'GGC', 'GCC', 'CCC']
        
        # make sure you convert to upper case first
        >>> reference = 'AaAtttGGgCcc'
        >>> mapper = ReadMapper(reference, 3)
        >>> mapper.sliding_window(reference, 3)
        ['AAA', 'AAT', 'ATT', 'TTT', 'TTG', 'TGG', 'GGG', 'GGC', 'GCC', 'CCC']
        
        """
        numSubSeq = (len(reference.upper()) - window_length)+1
        subSeqs = [reference.upper()[i:i+window_length] for i in range(0, numSubSeq)]
        newSubSeqs = [None if (re.search('N', x)) else x for x in subSeqs]
        return newSubSeqs

    def getAllIndices(self, string, sub, listIndices, offset):
        i = string.find(str(sub), offset)
        while i >= 0:
            listIndices.append(i)
            i = string.find(sub, i+1)
        return listIndices

    def build_index(self, reference, seed_length):
        """Create an index from all substrings in reference to their positions
        within reference. Returns a dictionary whose keys are the substrings of
        length window_length and whose values are a list of positions where that
        substring occurs.
        
        Note: you should use `sliding_window` to get all the substrings in `reference`.
        
        REMEMBER TO CONVERT REFERENCE TO UPPER CASE!
        
        For example:
        >>> build_index('AAATTTGGG', 2)
        {'AA':[0,1], 'AT':[2], 'TT':[3,4], 'TG':[5], 'GG':[6,7]}
        >>> build_index('AAATTTGGG', 3)
        {'AAA':[0], 'AAT':[1], 'ATT':[2], 'TTT':[3], 'TTG':[4], 'TGG':[5], 'GGG':[6]}
        
        """
        subStrings = self.sliding_window(reference.upper(), seed_length)
        subStringNoDup = []
        subStringNoDup = [i for i in subStrings if i not in subStringNoDup]
        indices = [self.getAllIndices(reference.upper(), x, [], 0) for x in subStringNoDup]
        self.ref_index = dict(zip(subStringNoDup, indices))
        for i in self.ref_index.keys():
            if i == None:
                del self.ref_index[i]
        return self.ref_index


    
