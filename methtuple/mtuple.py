from __future__ import print_function
import re
import itertools
import sys
import array

class MTuple:
    """ An MTuple instance stores methylation loci m-tuples and their associated counts for a single sample and a single "m" value.

    Attributes:
        sample_name: The name of the sample
        m: The "m" in m-tuples, i.e. size of the methylation-loci m-tuples.
        methylation_type: The type of methylation event: CG, CHG, CHH, CNN, CG/CHG, CG/CHH, CG/CNN, CHG/CHH, CHG/CNN, CHH/CNN, CG/CHG/CHH, CG/CHG/CNN, CG/CHH/CNN, CHG/CHH/CNN OR CG/CHG/CHH/CNN
        chr_map: A dictionary converting chromosome names to order, e.g. {'chr1': 1, 'chr10': 10, ..., 'chrX': 23}
        comethylation_patterns: A list of co-methylation patterns in alphabetic order, e.g. [MM, MU, UM, UU] for m = 2. Used as column names in output.
        counts: A dictionary storing positions and counts for each m-tuple. Keys are positions (stored as tuples) and values are their associated counts (stored as signed int arrays).
    """
    def __init__(self, sample_name, m, methylation_type, chr_map):
        """
        Initiates an MTuple instance with sample_name, m, methylation_type and chr_map given by arguments. The 'counts' field is set as an empty dictionary.
        """
        # Argument checks
        if not all([x in ['CG', 'CHG', 'CHH', 'CNN'] for x in methylation_type.split('/')]):
            raise ValueError("__init__ MTuple: 'methylation_type' must be one or more of 'CG', 'CHG', 'CHH' or 'CNN'. Multiple values must be separated by the '/' character.")
        if m < 0 or not isinstance(m, int):
            raise ValueError("__init__ MTuple: 'm' must be a positive integer.")
        # Initialise object of class MTuple
        self.sample_name = sample_name
        self.m = m
        self.methylation_type = methylation_type
        self.chr_map = chr_map
        tmp_k = list((itertools.product(('U', 'M'), repeat = m))) # Create all combinations of 'U', 'M' of length m by Cartesian product
        self.comethylation_patterns = sorted([''.join(a) for a in tmp_k])
        self.mtuples = {}
    def display(self):
        """Display an MTuple instance."""
        print('Sample name =', self.sample_name)
        print('Methylation type =', self.methylation_type)
        print('m = ', self.m)
        print('First 10 positions =', list(self.mtuples.keys())[:10])
        print('First 10 counts =', list(self.mtuples.values())[:10])
    def increment_count(self, pos, comethylation_state, read_1, read_2):
        """Increment the counts attribute based on the comethylation_state that has been extracted from read_1 and read_2.
        NB: read_2 should be set to None if data is single-end.
        NB: No check is made of the informative strand for read_1 and read_2."""
        # Collapse the 9 possible characters in the comethylation_state (., Z, z, X, x, H, h, U, u) to (., M, U, *unchanged*) based on the methylation_type parameter.
        if 'CG' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('Z', 'M')
            comethylation_state = comethylation_state.replace('z', 'U')
        if 'CHG' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('X', 'M')
            comethylation_state = comethylation_state.replace('x', 'U')
        if 'CHH' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('H', 'M')
            comethylation_state = comethylation_state.replace('h', 'U')
        if 'CNN' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('U', 'M')
            comethylation_state = comethylation_state.replace('u', 'U')

        # Check whether there are any unexpected, and therefore invalid, characters in comethylation_state
        if not re.search('[^MU]', comethylation_state) is None:
            exit_msg = ''.join([read_1.query_name, ' has an invalid comethylation state = ', comethylation_state, '.\nThis should never happen. Please file an issue at www.github.com/PeteHaitch/methtuple describing the error.'])
            sys.exit(exit_msg)
        if (len(comethylation_state) != self.m):
            exit_msg = ''.join(['Length of comethylation string (', str(len(comethylation_state)), ') does not equal m (', str(self.m), '). \nThis should never happen. Please file an issue at www.github.com/PeteHaitch/methtuple describing the error.'])
            sys.exit(exit_msg)

        # Check whether this m-tuple is already in MTuple. If it is just update the corresponding count, otherwise add that m-tuple and then update the corresponding count
        if not pos in self.mtuples:
            self.mtuples[pos] = array.array('i',[0] * 2**self.m)
        self.mtuples[pos][self.comethylation_patterns.index(comethylation_state)] += 1

__all__ = [
    'MTuple'
]
